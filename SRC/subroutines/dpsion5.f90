!> Calculate W-v zxqi(on the imaginary axis) and zxq(real axis) from sperctum weight rcxq.
! module m_getxc
!   implicit none
!   public :: getxc,mean_from_xc
!   private
! contains
!   ! Internal function: Compute the mean μ of a truncated normal distribution
!   ! given center xc and standard deviation sigma.
!   real(kind=8) function mean_from_xc(xc, sigma)
!     implicit none
!     real(kind=8), intent(in) :: xc, sigma
!     real(kind=8) :: Z, arg
!     ! Argument for the error function
!     arg = xc / (sqrt(2.0d0) * sigma)
!     ! Normalization constant over [0, ∞)
!     Z = sqrt(acos(-1.0d0)) * sigma / sqrt(2.0d0) * (1.0d0 + erf(arg))
!     ! Mean of the truncated normal distribution
!     mean_from_xc = xc + (sigma**2 / Z) * exp(-xc**2 / (2.0d0 * sigma**2))
!   end function mean_from_xc
!   ! Public function: Given mean μ and standard deviation σ,
!   ! numerically solve for the distribution center x_c.
!   real(kind=8) function getxc(mu, sigma)
!     use m_ftox
!     implicit none
!     real(kind=8), intent(in) :: mu, sigma
!     real(kind=8) :: a, b, c, fa, fb, fc, tol=1d-12
!     integer :: max_iter=100, i
!     ! Tolerance and iteration limit for bisection method
!     ! Initial bracket [a, b] for root finding
!     a =  max(0d0, -5.0d0 * sigma+mu)
!     b =  5.0d0 * sigma + mu
!     ! Evaluate residuals at endpoints
!     fa = mean_from_xc(a, sigma) - mu
!     fb = mean_from_xc(b, sigma) - mu
!     write (6,ftox) 'gggggggggg', a,b,sigma, 'xxxxxxx', ftod(mu),ftod(mean_from_xc(a, sigma)), ftod(mean_from_xc(b, sigma))
!     ! Check if root is bracketed
!     if (fa * fb > 0.0d0) then
!       getxc = -999.0d0  ! Return error value if no root found
!       return
!     end if
!     ! Bisection loop
!     do i = 1, max_iter
!       c = 0.5d0 * (a + b)
!       fc = mean_from_xc(c, sigma) - mu
!       if (abs(fc) < tol) then
!         getxc = c
!         return
!       end if
!       if (fa * fc < 0.0d0) then
!         b = c
!         fb = fc
!       else
!         a = c
!         fa = fc
!       end if
!     end do
!     ! Return approximate root if convergence not achieved
!     getxc = c
!   end function getxc
! end module m_getxc

module m_dpsion
  use m_kind, only: kp => kindrcxq
  use m_mpi,only: ipr
  public dpsion5, dpsion_init
  public dpsion_chiq_d, dpsion_chiq_h
  public dpsion_setup_rcxq_d, dpsion_setup_rcxq_h
  ! private
  real(8),allocatable :: his_L(:),his_R(:),his_C(:),rmat(:,:,:),rmatt(:,:,:),rmattx(:,:,:,:),imatt(:,:,:)
  complex(8),allocatable :: imattC(:,:,:)
  real(8),allocatable,save:: gfmat(:,:)
  logical:: eginit=.true., init=.true.
  complex(kind=kp), allocatable :: zxq_chipm(:,:,:)
#ifdef __GPU
  attributes (device) :: zxq_chipm
#endif
  complex(kind=kp), allocatable :: zxq_chipm_h(:,:,:)  ! host copy for dpsion_chiq_h
contains
  ! set omega-bin mesh, his_L, his_R, his_C, and Hilbert transformation weight rmat, rmatt, rmattx, imatt
  subroutine dpsion_init(realomega, imagomega, chipm)
    use m_freq, only:  frhis, freqr=>freq_r, freqi=>freq_i, nwhis, npm, nw_i, nw_w=>nw, niwt=>niw
    use m_lgunit,only:stdo
    implicit none
    logical, intent(in)::realomega, imagomega, chipm
    complex(8):: img=(0d0,1d0), zz, rrr(-nwhis:nwhis)
    real(8),parameter:: pi  = 4d0*datan(1d0)
    integer :: it
    if(.not.init) return
    allocate( his_L(-nwhis:nwhis),source=[-frhis(nwhis+1:1+1:-1),0d0,frhis(1  :nwhis)  ])
    allocate( his_R(-nwhis:nwhis),source=[-frhis(nwhis  :1  :-1),0d0,frhis(1+1:nwhis+1)])
    allocate( his_C(-nwhis:nwhis),source=(his_L+his_R)/2d0) !bins are [his_Left,his_Right] !his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
    realomegacase: if(realomega)then
      if(ipr) write(stdo,*) " --- realomega --- "
      if(npm==1) then
        allocate(rmat(0:nw_w,-nwhis:nwhis,npm), source=0d0)
        do it =  0, nw_w
          zz = freqr(it)
          call hilbertmat(zz,  nwhis,his_L,his_C,his_R, rrr)
          rmat(it,:,1) = dreal(rrr)/pi
        enddo
        if(chipm) then
          allocate( rmattx(0:nw_w,nwhis,npm,2) )
          rmattx(:,1:nwhis,1,1) =  rmat(:,1:nwhis,1)
          rmattx(:,1:nwhis,1,2) = -rmat(:,-1:-nwhis:-1,1)
        else  
          allocate( rmatt(0:nw_w,nwhis,npm) )
          rmatt(:,1:nwhis,1) =  rmat(:,1:nwhis,1) - rmat(:,-1:-nwhis:-1,1)
        endif
        deallocate(rmat)
      elseif(npm==2) then
        allocate(rmatt(-nw_w:nw_w,nwhis,npm))
        do it  =  -nw_w,nw_w
          zz = merge(-freqr(-it),freqr(it),it<0) 
          call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
          rmatt(it,:,1) =  dreal(rrr ( 1: nwhis))/pi
          rmatt(it,:,2) = -dreal(rrr(-1:-nwhis:-1))/pi
        enddo
      endif
    endif realomegacase
    imagomecacase: if(imagomega) then
      write(stdo,*) " --- imagomega --- "
      if(npm==1) then
        allocate( imatt(niwt, nwhis,npm) )
        do it =  1,niwt
          zz = img*freqi(it)  
          call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
          imatt(it,1:nwhis,1) = dreal(rrr(1:nwhis) - rrr(-1:-nwhis:-1))/pi
        enddo
      else ! npm=2 case 
        allocate( imattC(niwt, nwhis,npm) )
        do it =  1,niwt
          zz = img*freqi(it)  
          call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
          imattC(it,1:nwhis,1) =   rrr( 1: nwhis   )/pi
          imattC(it,1:nwhis,2) = - rrr(-1:-nwhis:-1)/pi
        enddo
      endif
    endif imagomecacase
    init=.false.
  end subroutine dpsion_init
  subroutine dpsion_chiq_d(realomega, imagomega, chipm, rcxq, zxqi, npr, npr_col, schi, isp, ecut)
    use m_keyvalue,only: getkeyvalue
    use m_GWinput, only: gwinput_init, gwinput_loaded, tg_SmearX0 => SmearX0
    use m_freq, only: frhis, freqr=>freq_r,freqi=>freq_i, nwhis, npm, nw_i, nw_w=>nw, niwt=>niw
!    use m_readgwinput, only: egauss
    use m_ftox
    use m_lgunit, only: stdo
    use m_blas, only: m_op_T
#if defined(__MP) && defined(__GPU)
    use m_blas, only: gemm => cmm_d
#elif defined(__MP)
    use m_blas, only: gemm => cmm_h
#elif defined(__GPU)
    use m_blas, only: gemm => zmm_d
#else
    use m_blas, only: gemm => zmm_h
#endif
    implicit none
    logical, intent(in):: realomega, imagomega, chipm
    real(8), intent(in):: ecut, schi
    integer, intent(in):: isp, npr, npr_col
    complex(kind=kp), intent(inout):: rcxq(1:npr,1:npr_col,(1-npm)*nwhis:nwhis)
    complex(kind=kp), intent(out) ::  zxqi(1:npr,1:npr_col,niwt)
    complex(kind=kp), parameter:: CONE = (1_kp, 0_kp), CZERO = (0_kp, 0_kp)
    integer :: iw,i,j
    real(8), parameter:: pi  = 4d0*datan(1d0)
    complex(8), parameter :: img = (0d0,1d0)
    complex(kind=kp) :: zxq_work(1:npr,nw_i:nw_w), cimatt(niwt,nwhis,npm), crmatt(nw_i:nw_w,nwhis,npm)
    complex(kind=kp), allocatable :: rcxq_work(:,:), cgfmat(:,:)
    integer :: ipr_col, ipm, istat, ispx
    real(8) :: wfac,smearx0
#ifdef __GPU
    attributes(device) :: rcxq, zxqi
#endif
    if(ipr) write(stdo,ftox)" -- dpsion_chiq_d:start... nw_w nwhis=",nw_w,nwhis
    call flush(stdo)
    if(chipm.and.npm==2) call rx( 'x0kf_v4h:npm==2 .AND. chipm is not meaningful probably')  ! Note rcxq here is negative 
    !$acc data copyin(his_R, his_L)
    call gwinput_init()
    if (gwinput_loaded) then
       smearx0 = tg_SmearX0
    else
       call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
!       call getkeyvalue("GWinput","SmearX0", smearx0, default=0d0 )
    endif
    GaussianFilter: if(abs(smearx0)>1d-15) then
      if(ipr) write(6,'("SmearX0= ",d13.6)') smearx0
      allocate(gfmat(nwhis,nwhis))
      allocate(cgfmat(nwhis,nwhis))
      allocate(rcxq_work(npr,nwhis))
      if(ipr) write(stdo,ftox) 'dpsion_chiq: SmearX0 (chi0 GaussianFilter) is not checked yet: see dpsion_chiq'
      gfmat=gaussianfilterhis(smearx0,frhis,nwhis)

      !$acc data copyin(gfmat) create(cgfmat, rcxq_work)
      !$acc kernels
      cgfmat(:,:) = cmplx(gfmat(:,:), kind=kp)
      !$acc end kernels
      do ipr_col = 1, npr_col
        !$acc kernels
        rcxq_work(1:npr,1:nwhis) = rcxq(1:npr,ipr_col,1:nwhis)
        !$acc end kernels
        istat = gemm(rcxq_work, cgfmat, rcxq(1,ipr_col,1), npr, nwhis, nwhis, ldC=npr*npr_col, opB=m_op_T)
      enddo
      if(npm==2) then
        !$acc kernels
        cgfmat(1:nwhis,1:nwhis) = cmplx(gfmat(1:nwhis,nwhis:1:-1), kind=kp)
        !$acc end kernels
        do ipr_col = 1, npr_col
          !$acc kernels
          rcxq_work(1:npr,1:nwhis) = rcxq(1:npr,ipr_col,-nwhis:-1:1)
          !$acc end kernels
          istat = gemm(rcxq_work, cgfmat, rcxq(1,ipr_col,-nwhis), npr, nwhis, nwhis, ldC=npr*npr_col, opB=m_op_T)
        enddo
      endif
      !$acc end data
      deallocate(gfmat)
    endif GaussianFilter
      
    ispx = merge(isp,3-isp,schi>=0) !  if(schi<0)  ispx = 3-isp  
    if(realomega.and.nwhis <= nw_w) call rxii('dpsion5: nwhis<=nw_w',nwhis,nw_w)
    if(realomega.and.freqr(0)/=0d0) call rx( 'dpsion5: freqr(0)/=0d0') ! I think current version allows any freqr(iw), independent from frhis.

    call flush(stdo)
    do iw= 1, nwhis
      wfac=merge(exp(-(his_C(iw)/ecut)**2 ),1d0, ecut<1d9)     ! rcxq= Average value of Im chi.    Note rcxq is "negative" (
      !$acc kernels
      rcxq(1:npr,1:npr_col,iw)= -wfac/(his_R(iw)-his_L(iw))*rcxq(1:npr,1:npr_col,iw)
      !$acc end kernels
    enddo
    if(npm==2) then ! 2025-12-05 Bugfix for npm=2 case
      do iw= -nwhis, -1
        wfac=merge(exp(-(his_C(iw)/ecut)**2 ),1d0, ecut<1d9)     ! rcxq= Average value of Im chi.    Note rcxq is "negative" (
        !$acc kernels
        rcxq(1:npr,1:npr_col,iw)= -wfac/(his_R(iw)-his_L(iw))*rcxq(1:npr,1:npr_col,iw)
        !$acc end kernels
      enddo
    endif
    !$acc end data
    if_IMAGOMEGA: if(imagomega) then !Hilbert Transformation to get real part
      if(ipr) write(stdo,ftox)" -- dpsion_chiq_d:start imagomega"
      if(npm==1) then
        !$acc data copyin(imatt) create(cimatt)
        !$acc kernels
        cimatt(:,:,:) = cmplx(imatt(:,:,:), kind=kp)
        !$acc end kernels
        ! istat = gemm(rcxq(1,1,1), cimatt, zxqi, npr*npr_col, niwt, nwhis, opB=m_op_T)
        ! Above line is replaced by the following loop to reduce internal memory usage on gemmul8
        do ipr_col = 1, npr_col
          istat = gemm(rcxq(1,ipr_col,1), cimatt, zxqi(1,ipr_col,1), m=npr, n=niwt, k=nwhis, &
                     & ldA=npr*npr_col, opB=m_op_T, ldC=npr*npr_col)
        enddo
        !$acc end data
      elseif(npm==2) then
        !$acc data copyin(imattC) create(cimatt)
        !$acc kernels
        cimatt(:,1:nwhis,1) = cmplx(imattC(:,1:nwhis: 1,1), kind=kp)
        cimatt(:,1:nwhis,2) = cmplx(imattC(:,nwhis:1:-1,2), kind=kp)
        !$acc end kernels
        istat = gemm(rcxq(1,1,     1), cimatt(:,1,1), zxqi, npr*npr_col, niwt, nwhis, opB=m_op_T)
        istat = gemm(rcxq(1,1,-nwhis), cimatt(:,1,2), zxqi, npr*npr_col, niwt, nwhis, opB=m_op_T, beta=CONE)
        !$acc end data
      endif
      if(ipr) write(stdo,ftox)" -- dpsion_chiq_d:end of imagomega"
    endif if_IMAGOMEGA
    if_REALOMEGA: if(realomega) then !Hilbert Transformation to get real part
      if(ipr) write(stdo,ftox)" -- dpsion_chiq_d:start realomega"
      if(npm == 1 .and. .not.chipm) then
        !$acc data copyin(rmatt) create(crmatt, zxq_work)
        !$acc kernels
        crmatt(:,:,:) = cmplx(rmatt(:,:,:), kind=kp)
        !$acc end kernels
        do ipr_col = 1, npr_col
          istat = gemm(rcxq(1,ipr_col,1), crmatt, zxq_work, npr, nw_w+1, nwhis, ldA=npr*npr_col, opB=m_op_T)
          !$acc kernels
          rcxq(1:npr,ipr_col,0:nw_w) = rcxq(1:npr,ipr_col,0:nw_w)*img + zxq_work(1:npr,0:nw_w)
          !$acc end kernels
        enddo
        !$acc end data
      elseif(npm == 1 .and. chipm) then
        if(.not.allocated(zxq_chipm)) then
          allocate(zxq_chipm(npr,npr_col,nw_i:nw_w))
          !$acc kernels
          zxq_chipm(:,:,:) = (0_kp, 0_kp)
          !$acc end kernels
        endif
        if(ispx == 1) then
          !$acc kernels
          zxq_chipm(:,:,1:nw_w) = zxq_chipm(:,:,1:nw_w)+ img*rcxq(:,:,1:nw_w) 
          !$acc end kernels
        endif
        !$acc data copyin(rmattx) create(crmatt)
        !$acc kernels
        crmatt(:,:,:) = cmplx(rmattx(:,:,:,ispx), kind=kp)
        !$acc end kernels
        istat = gemm(rcxq(1,1,1), crmatt, zxq_chipm, npr*npr_col, nw_w+1, nwhis, opB=m_op_T, beta=CONE)
        !$acc end data
      elseif(npm == 2) then
        !$acc data copyin(rmatt) create(crmatt, zxq_work)
        !$acc kernels
        crmatt(:,1:nwhis,1) = cmplx(rmatt(:,1:nwhis: 1,1), kind=kp)
        crmatt(:,1:nwhis,2) = cmplx(rmatt(:,nwhis:1:-1,2), kind=kp)
        !$acc end kernels
        do ipr_col = 1, npr_col
          istat = gemm(rcxq(1,ipr_col,     1), crmatt(:,1,1), zxq_work, npr, (nw_w-nw_i)+1, nwhis, ldA=npr*npr_col, opB=m_op_T)
          istat = gemm(rcxq(1,ipr_col,-nwhis), crmatt(:,1,2), zxq_work, npr, (nw_w-nw_i)+1, nwhis, ldA=npr*npr_col, opB=m_op_T,&
                       beta=CONE)
          !$acc kernels
          rcxq(1:npr,ipr_col,nw_i:nw_w) = rcxq(1:npr,ipr_col,nw_i:nw_w)*img + zxq_work(1:npr,nw_i:nw_w) !override
          !$acc end kernels
        enddo
        !$acc end data
      endif
      if(ipr) write(stdo,ftox)" -- dpsion_chiq_d:end of realomega"
    endif if_REALOMEGA
    call flush(stdo)
  end subroutine dpsion_chiq_d

  subroutine dpsion_setup_rcxq_d(rcxq, npr, npr_col, isp)
    use m_freq, only:nwhis, npm, nw_i, nw_w => nw
    implicit none
    integer, intent(in) :: npr, npr_col, isp
    complex(kind=kp), intent(inout):: rcxq(1:npr,1:npr_col,(1-npm)*nwhis:nwhis)
    if(isp == 1) then
      !$acc kernels
      rcxq(:,:,:) = (0_kp, 0_kp)
      !$acc end kernels
    elseif(isp == 2) then
      if(allocated(zxq_chipm)) then
        !$acc kernels
        rcxq(:,:,nw_i:nw_w) = zxq_chipm(:,:,nw_i:nw_w)
        !$acc end kernels
        deallocate(zxq_chipm)
      endif
    endif
  end subroutine dpsion_setup_rcxq_d

  !> Host (CPU) version of dpsion_setup_rcxq: uses zxq_chipm_h, no OpenACC.
  subroutine dpsion_setup_rcxq_h(rcxq, npr, npr_col, isp)
    use m_freq, only: nwhis, npm, nw_i, nw_w => nw
    implicit none
    integer, intent(in) :: npr, npr_col, isp
    complex(kind=kp), intent(inout) :: rcxq(1:npr,1:npr_col,(1-npm)*nwhis:nwhis)
    if (isp == 1) then
      rcxq(:,:,:) = (0_kp, 0_kp)
    elseif (isp == 2) then
      if (allocated(zxq_chipm_h)) then
        rcxq(:,:,nw_i:nw_w) = zxq_chipm_h(:,:,nw_i:nw_w)
        deallocate(zxq_chipm_h)
      endif
    endif
  end subroutine dpsion_setup_rcxq_h

  !> Host (CPU) version of dpsion_chiq: no attributes(device), no OpenACC.
  !> rcxq and zxqi must be host-resident (e.g. pointing to shm_wvr/shm_wvi).
  subroutine dpsion_chiq_h(realomega, imagomega, chipm, rcxq, zxqi, npr, npr_col, schi, isp, ecut, smearx0_in)
    use m_keyvalue, only: getkeyvalue
    use m_GWinput, only: gwinput_init, gwinput_loaded, tg_SmearX0 => SmearX0
    use m_freq, only: frhis, freqr=>freq_r, freqi=>freq_i, nwhis, npm, nw_i, nw_w=>nw, niwt=>niw
    use m_ftox
    use m_lgunit, only: stdo
    use m_blas, only: m_op_T
#if defined(__MP)
    use m_blas, only: gemm => cmm_h
#else
    use m_blas, only: gemm => zmm_h
#endif
    implicit none
    logical, intent(in)  :: realomega, imagomega, chipm
    real(8), intent(in)  :: ecut, schi
    real(8), intent(in), optional :: smearx0_in  ! caller-supplied SmearX0 (e.g. SmearX0q0 at offset-Gamma); falls back to toml SmearX0 if absent
    integer, intent(in)  :: isp, npr, npr_col
    complex(kind=kp), intent(inout) :: rcxq(1:npr,1:npr_col,(1-npm)*nwhis:nwhis)
    complex(kind=kp), intent(out)   :: zxqi(1:npr,1:npr_col,niwt)
    complex(kind=kp), parameter :: CONE = (1_kp, 0_kp), CZERO = (0_kp, 0_kp)
    integer  :: iw, i, j
    real(8),  parameter :: pi = 4d0*datan(1d0)
    complex(8), parameter :: img = (0d0, 1d0)
    complex(kind=kp) :: zxq_work(1:npr,nw_i:nw_w), cimatt(niwt,nwhis,npm), crmatt(nw_i:nw_w,nwhis,npm)
    complex(kind=kp), allocatable :: rcxq_work(:,:), cgfmat(:,:)
    integer  :: ipr_col, ipm, istat, ispx
    real(8)  :: wfac, smearx0
    real(8)  :: fs_b(2), ar_b(2), fs_a(2), ar_a(2), frcw   ! SmearX0 f-sum/area diagnostic (head & off-diag)
    if (ipr) write(stdo,ftox) " -- dpsion_chiq_h: start... nw_w nwhis=", nw_w, nwhis
    call flush(stdo)
    if (chipm.and.npm==2) call rx('dpsion_chiq_h: npm==2 .AND. chipm is not meaningful')
    call gwinput_init()
    if (present(smearx0_in)) then
      smearx0 = smearx0_in
    elseif (gwinput_loaded) then
      smearx0 = tg_SmearX0
    else
      call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
    endif
    GaussianFilter: if (abs(smearx0) > 1d-15) then
      if (ipr) write(6,'("SmearX0(chi0 GaussianFilter)= ",d13.6)') smearx0
      allocate(gfmat(nwhis,nwhis))
      allocate(cgfmat(nwhis,nwhis))
      allocate(rcxq_work(npr,nwhis))
      gfmat = gaussianfilterhis(smearx0, frhis, nwhis)
      cgfmat(:,:) = cmplx(gfmat(:,:), kind=kp)
      ! --- f-sum/area diagnostic BEFORE filtering (head=(1,1), off-diag=(1,2)) ---
      ar_b=0d0; fs_b=0d0
      do iw=1,nwhis
        frcw=(frhis(iw)+frhis(iw+1))/2d0
        ar_b(1)=ar_b(1)+dble(rcxq(1,1,iw)); fs_b(1)=fs_b(1)+frcw*dble(rcxq(1,1,iw))
        if(npr_col>=2) then
          ar_b(2)=ar_b(2)+dble(rcxq(1,2,iw)); fs_b(2)=fs_b(2)+frcw*dble(rcxq(1,2,iw))
        endif
      enddo
      do ipr_col = 1, npr_col
        rcxq_work(1:npr,1:nwhis) = rcxq(1:npr,ipr_col,1:nwhis)
        istat = gemm(rcxq_work, cgfmat, rcxq(1,ipr_col,1), npr, nwhis, nwhis, ldC=npr*npr_col, opB=m_op_T)
      enddo
      ! --- f-sum/area diagnostic AFTER filtering; f-sum must be (nearly) unchanged ---
      ar_a=0d0; fs_a=0d0
      do iw=1,nwhis
        frcw=(frhis(iw)+frhis(iw+1))/2d0
        ar_a(1)=ar_a(1)+dble(rcxq(1,1,iw)); fs_a(1)=fs_a(1)+frcw*dble(rcxq(1,1,iw))
        if(npr_col>=2) then
          ar_a(2)=ar_a(2)+dble(rcxq(1,2,iw)); fs_a(2)=fs_a(2)+frcw*dble(rcxq(1,2,iw))
        endif
      enddo
      write(stdo,"(' SmearX0 fsum-diag head(1,1):  area b/a=',2es13.5,'  f-sum b/a=',2es13.5)") ar_b(1),ar_a(1),fs_b(1),fs_a(1)
      if(npr_col>=2) write(stdo,"(' SmearX0 fsum-diag offd(1,2): area b/a=',2es13.5,'  f-sum b/a=',2es13.5)") ar_b(2),ar_a(2),fs_b(2),fs_a(2)
      if (npm == 2) then
        cgfmat(1:nwhis,1:nwhis) = cmplx(gfmat(1:nwhis,nwhis:1:-1), kind=kp)
        do ipr_col = 1, npr_col
          rcxq_work(1:npr,1:nwhis) = rcxq(1:npr,ipr_col,-nwhis:-1:1)
          istat = gemm(rcxq_work, cgfmat, rcxq(1,ipr_col,-nwhis), npr, nwhis, nwhis, ldC=npr*npr_col, opB=m_op_T)
        enddo
      endif
      deallocate(gfmat, cgfmat, rcxq_work)
    endif GaussianFilter

    ispx = merge(isp, 3-isp, schi >= 0)
    if (realomega .and. nwhis <= nw_w) call rxii('dpsion_chiq_h: nwhis<=nw_w', nwhis, nw_w)
    if (realomega .and. freqr(0)/=0d0) call rx('dpsion_chiq_h: freqr(0)/=0d0')
    call flush(stdo)
    do iw = 1, nwhis
      wfac = merge(exp(-(his_C(iw)/ecut)**2), 1d0, ecut < 1d9)
      rcxq(1:npr,1:npr_col,iw) = -wfac/(his_R(iw)-his_L(iw)) * rcxq(1:npr,1:npr_col,iw)
    enddo
    if (npm == 2) then
      do iw = -nwhis, -1
        wfac = merge(exp(-(his_C(iw)/ecut)**2), 1d0, ecut < 1d9)
        rcxq(1:npr,1:npr_col,iw) = -wfac/(his_R(iw)-his_L(iw)) * rcxq(1:npr,1:npr_col,iw)
      enddo
    endif

    if_IMAGOMEGA: if (imagomega) then
      if (ipr) write(stdo,ftox) " -- dpsion_chiq_h: start imagomega"
      if (npm == 1) then
        cimatt(:,:,:) = cmplx(imatt(:,:,:), kind=kp)
        do ipr_col = 1, npr_col
          istat = gemm(rcxq(1,ipr_col,1), cimatt, zxqi(1,ipr_col,1), m=npr, n=niwt, k=nwhis, &
                     & ldA=npr*npr_col, opB=m_op_T, ldC=npr*npr_col)
        enddo
      elseif (npm == 2) then
        cimatt(:,1:nwhis,1) = cmplx(imattC(:,1:nwhis: 1,1), kind=kp)
        cimatt(:,1:nwhis,2) = cmplx(imattC(:,nwhis:1:-1,2), kind=kp)
        istat = gemm(rcxq(1,1,     1), cimatt(:,1,1), zxqi, npr*npr_col, niwt, nwhis, opB=m_op_T)
        istat = gemm(rcxq(1,1,-nwhis), cimatt(:,1,2), zxqi, npr*npr_col, niwt, nwhis, opB=m_op_T, beta=CONE)
      endif
      if (ipr) write(stdo,ftox) " -- dpsion_chiq_h: end of imagomega"
    endif if_IMAGOMEGA

    if_REALOMEGA: if (realomega) then
      if (ipr) write(stdo,ftox) " -- dpsion_chiq_h: start realomega"
      if (npm == 1 .and. .not.chipm) then
        crmatt(:,:,:) = cmplx(rmatt(:,:,:), kind=kp)
        do ipr_col = 1, npr_col
          istat = gemm(rcxq(1,ipr_col,1), crmatt, zxq_work, npr, nw_w+1, nwhis, ldA=npr*npr_col, opB=m_op_T)
          rcxq(1:npr,ipr_col,0:nw_w) = rcxq(1:npr,ipr_col,0:nw_w)*img + zxq_work(1:npr,0:nw_w)
        enddo
      elseif (npm == 1 .and. chipm) then
        if (.not.allocated(zxq_chipm_h)) then
          allocate(zxq_chipm_h(npr,npr_col,nw_i:nw_w))
          zxq_chipm_h(:,:,:) = (0_kp, 0_kp)
        endif
        if (ispx == 1) then
          zxq_chipm_h(:,:,1:nw_w) = zxq_chipm_h(:,:,1:nw_w) + img*rcxq(:,:,1:nw_w)
        endif
        crmatt(:,:,:) = cmplx(rmattx(:,:,:,ispx), kind=kp)
        istat = gemm(rcxq(1,1,1), crmatt, zxq_chipm_h, npr*npr_col, nw_w+1, nwhis, opB=m_op_T, beta=CONE)
      elseif (npm == 2) then
        crmatt(:,1:nwhis,1) = cmplx(rmatt(:,1:nwhis: 1,1), kind=kp)
        crmatt(:,1:nwhis,2) = cmplx(rmatt(:,nwhis:1:-1,2), kind=kp)
        do ipr_col = 1, npr_col
          istat = gemm(rcxq(1,ipr_col,     1), crmatt(:,1,1), zxq_work, npr, (nw_w-nw_i)+1, nwhis, ldA=npr*npr_col, opB=m_op_T)
          istat = gemm(rcxq(1,ipr_col,-nwhis), crmatt(:,1,2), zxq_work, npr, (nw_w-nw_i)+1, nwhis, ldA=npr*npr_col, opB=m_op_T, beta=CONE)
          rcxq(1:npr,ipr_col,nw_i:nw_w) = rcxq(1:npr,ipr_col,nw_i:nw_w)*img + zxq_work(1:npr,nw_i:nw_w)
        enddo
      endif
      if (ipr) write(stdo,ftox) " -- dpsion_chiq_h: end of realomega"
    endif if_REALOMEGA
    call flush(stdo)
  end subroutine dpsion_chiq_h

  subroutine dpsion5(realomega,imagomega,rcxq,nmbas1,nmbas2, zxq,zxqi, chipm,schi,isp,ecut,ecuts)
    use m_freq,only:  frhis, freqr=>freq_r,freqi=>freq_i, nwhis, npm, nw_i, nw_w=>nw, niwt=>niw
!    use m_readgwinput,only: egauss
!    use m_GaussianFilter,only: GaussianFilter
    use m_ftox
    use m_lgunit,only:stdo
    use m_kind,only:kindrcxq
    use m_keyvalue,only: getkeyvalue
    use m_GWinput, only: gwinput_init, gwinput_loaded, tg_SmearX0 => SmearX0
    implicit none
    intent(in)::     realomega,imagomega,     nmbas1,nmbas2,           chipm,schi,isp,ecut,ecuts
    intent(out)::                        rcxq,                zxq,zxqi
    !                                    rcxq is destroyed
    !  works for timereversal=F (npm=2 case).
    !input
    !i   frhis(1:nwhis+1) : specify histgram bins i-th bin is [frhis(i), frhis(i+1)].
    !i   rcxq: the spectrum weight for given bins along the real-axis.
    !i   freqr (0:nw_w) : Calcualte zxq for these real energies.
    !i   freqi (1:niwt) : Calcualte zxqi for these imaginary energies.
    !i   realomega  : A switch to calculate zxq or not.
    !i   imagomega: : A switch to calculate zxqi or not.
    !o   zxq:  W-v along the real axis on freqr(0:nw_w). not accumlating
    !o   zxqi: W-v along the imag axis on freqi(niwt). not accumlating
    !r  We suppose "freqr(i)=moddle of i-th bin; freqr(0)=0." (I think called routine hilbertmat itself is not limited by this condition).
    integer:: igb1,igb2, iw,iwp,ix,ifxx,nmbas1,nmbas2,isp,ispx,it, ii,i,ibas1,ibas2,nmnm
    logical :: evaltest     
    real(8):: px,omp,om,om2,om1, aaa,d_omg, ecut,ecuts,wcut,dee,schi, domega_r,domega_c,domega_l,delta_l,delta_r,smearx0
    complex(8):: zxq(nmbas1,nmbas2, nw_i:nw_w),zxqi(nmbas1,nmbas2,niwt),img=(0d0,1d0),beta,wfac, zz,rrr(-nwhis:nwhis)
    logical :: realomega, imagomega,chipm,debug=.false.
    integer:: jpm,ipm,verbose,isgi   !     complex(8):: x0mean(nw_i:nw_w,nmbas,nmbas)
    real(8),parameter:: pi  = 4d0*datan(1d0)
    logical::init=.true.
    integer:: imbas1,imbas2,j
    complex(8):: rcxqin(1:nwhis)
    complex(kindrcxq):: rcxq(nmbas1,nmbas2, nwhis,npm)

    if(ipr) write(stdo,ftox)" -- dpsion5: start... nw_w nwhis=",nw_w,nwhis
    if(chipm.and.npm==2) call rx( 'x0kf_v4h:npm==2 .AND. chipm is not meaningful probably')  ! Note rcxq here is negative 
    call cputid(0)
    call gwinput_init()
    if (gwinput_loaded) then
       smearx0 = tg_SmearX0
    else
       call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
!       call getkeyvalue("GWinput","SmearX0", smearx0, default=0d0 )
    endif
    GaussianFilter: if(abs(smearx0)>1d-15) then
       if(eginit) then
         if(ipr) write(stdo,'("SmearX0= ",d13.6)') smearx0
          allocate(gfmat(nwhis,nwhis))
          gfmat=gaussianfilterhis(smearx0,frhis,nwhis)
          eginit=.false.
       endif
       do ipm=1,npm
          do imbas1=1,nmbas1
             do imbas2=1,nmbas2
                rcxqin = rcxq(imbas1,imbas2,1:nwhis,ipm)
                rcxq(imbas1,imbas2,1:nwhis,ipm) = matmul(gfmat,rcxqin)
             enddo
          enddo
       enddo       !write(6,"(' End of Gaussian Filter egauss=',f9.4)") egauss
    endif GaussianFilter
       
    ispx = merge(isp,3-isp,schi>=0) !  if(schi<0)  ispx = 3-isp  
    if(realomega.and.nwhis <= nw_w) call rxii('dpsion5: nwhis<=nw_w',nwhis,nw_w)
    if(realomega.and.freqr(0)/=0d0) call rx( 'dpsion5: freqr(0)/=0d0') ! I think current version allows any freqr(iw), independent from frhis.

    if(init) then !get Hilbert transformation matrix
      allocate( his_L(-nwhis:nwhis),source=[-frhis(nwhis+1:1+1:-1),0d0,frhis(1  :nwhis)  ])
      allocate( his_R(-nwhis:nwhis),source=[-frhis(nwhis  :1  :-1),0d0,frhis(1+1:nwhis+1)])
      allocate( his_C(-nwhis:nwhis),source=(his_L+his_R)/2d0) !bins are [his_Left,his_Right] !his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
      realomegacase: if(realomega)then;     if(ipr) write(stdo,*) " --- realomega --- "
        if(npm==1) then
          allocate(rmat(0:nw_w,-nwhis:nwhis,npm), source=0d0)
          do it =  0,nw_w
            zz = freqr(it)
            call hilbertmat(zz,  nwhis,his_L,his_C,his_R, rrr)
            rmat(it,:,1) = dreal(rrr)/pi
          enddo
          if(chipm) then
            allocate( rmattx(0:nw_w,nwhis,npm,2) )
            rmattx(:,1:nwhis,1,1) =  rmat(:,1:nwhis,1)
            rmattx(:,1:nwhis,1,2) = -rmat(:,-1:-nwhis:-1,1)
          else  
            allocate( rmatt(0:nw_w,nwhis,npm) )
            rmatt(:,1:nwhis,1) =  rmat(:,1:nwhis,1) - rmat(:,-1:-nwhis:-1,1)
          endif
          deallocate(rmat)
        else  ! npm==2 
          allocate(rmatt(-nw_w:nw_w,nwhis,npm))
          do it  =  -nw_w,nw_w
            zz = merge(-freqr(-it),freqr(it),it<0) 
            call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
            rmatt(it,:,1) =  dreal(rrr  (1:nwhis))/pi
            rmatt(it,:,2) = -dreal(rrr(-1:-nwhis:-1))/pi
          enddo
        endif
      endif realomegacase
      imagomecacase: if(imagomega) then
        if(npm==1) then
          allocate( imatt(niwt, nwhis,npm) )
          do it =  1,niwt
            zz = img*freqi(it)  
            call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            imatt(it,1:nwhis,1) = dreal(rrr(1:nwhis) - rrr(-1:-nwhis:-1))/pi
          enddo
        else ! npm=2 case 
          allocate( imattC(niwt, nwhis,npm) )
          do it =  1,niwt
            zz = img*freqi(it)  
            call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            imattC(it,1:nwhis,1) =   rrr( 1: nwhis   )/pi
            imattC(it,1:nwhis,2) = - rrr(-1:-nwhis:-1)/pi
          enddo
        endif
      endif imagomecacase
      init=.false.
    endif

    do iw= 1, nwhis
      wfac=merge(exp(-(his_C(iw)/ecut)**2 ),1d0, ecut<1d9)     ! rcxq= Average value of Im chi.    Note rcxq is "negative" (
      rcxq(:,:,iw,:)= -wfac/(his_r(iw)-his_l(iw))*rcxq(:,:,iw,:)
    enddo
    if(realomega) then !Hilbert Transformation to get real part
      if(chipm.AND.ispx==1) zxq(:,:,1:nw_w)= zxq(:,:,1:nw_w)+ img*rcxq(:,:,1:nw_w,1) 
      if(.not.chipm)        zxq(:,:,1:nw_w)= img*rcxq(:,:,1:nw_w,1)
      nmnm=2*nmbas1*nmbas2
      if(npm==1.and.chipm) then
        call dgemm('n','t',nmnm,nw_w+1,    nwhis,1d0,dcmplx(rcxq),         nmnm,rmattx(:,:,:,ispx),nw_w+1,1d0,zxq,nmnm)
      elseif(npm==1) then
        call dgemm('n','t',nmnm,nw_w+1,    nwhis,1d0,dcmplx(rcxq),         nmnm,rmatt,             nw_w+1,1d0,zxq,nmnm)
      elseif(npm==2) then
         zxq(:,:,-1:-nw_w:-1)=zxq(:,:,-1:-nw_w:-1) + img*rcxq(:,:,1:nw_w,2)
         !call zaxpy( nmbas1*nmbas2, img, rcxq(1,1,iw,2),1, zxq(:,:,-iw),1)
        call dgemm('n','t',nmnm,npm*nw_w+1,nwhis,1d0,dcmplx(rcxq(:,:,:,1)),nmnm,rmatt(:,:,1),npm*nw_w+1,1d0,zxq,nmnm)
        call dgemm('n','t',nmnm,npm*nw_w+1,nwhis,1d0,dcmplx(rcxq(:,:,:,2)),nmnm,rmatt(:,:,2),npm*nw_w+1,1d0,zxq,nmnm)
      endif
    endif
    if(imagomega) then !Hilbert Transformation to get real part
      nmnm=nmbas1*nmbas2
      if(npm==1) then
        call dgemm('n','t',2*nmnm,niwt,nwhis,1d0,dcmplx(rcxq), 2*nmnm, imatt, niwt, 0d0, zxqi, 2*nmnm )
      elseif(npm==2) then
        call zgemm('n','t', nmnm,niwt,nwhis,1d0, dcmplx(rcxq(:,:,:,1)),nmnm,imattC(1,1,1),niwt, 0d0,zxqi, nmnm )
        call zgemm('n','t', nmnm,niwt,nwhis,1d0, dcmplx(rcxq(:,:,:,2)),nmnm,imattC(1,1,2),niwt, 1d0,zxqi, nmnm )
      endif
    endif
    if(ipr) write(stdo,'("         end dpsion5 ",$)')
    call cputid(0)
  end subroutine dpsion5
!  subroutine GaussianFilter(rcxq,nmbas1,nmbas2, egauss,iprint)
  function gaussianfilterhis(smearx0, frhis,nwhis) result(gfmat)
    ! Bin-width-weighted, weight-conserving Gaussian smoothing of Im chi0 along omega.
    ! frhis is a strongly non-uniform (exponential) mesh, so the kernel MUST carry the bin
    ! width dfr; a plain count-sum normalization (the old code) biased toward the dense
    ! omega~0 bins and broke the sum rules.
    ! gemm computes rcxq_new(i) = sum_j rcxq(j) * gfmat(i,j)  (target i, source j).
    ! We use the scatter form  gfmat(i,j) = K(i,j) dfr(i) / sum_i' K(i',j) dfr(i')  so that
    !   sum_i gfmat(i,j) = 1  for every source j
    ! => area  sum_iw rcxq(iw)        is conserved EXACTLY (0th moment), and
    !    f-sum sum_iw omega_iw rcxq   is conserved to the continuum (symmetric-kernel) limit.
    ! Reduces to the identity as smearx0 -> 0.
    implicit none
    integer,intent(in):: nwhis
    real(8),intent(in):: smearx0,frhis(nwhis+1)
    real(8):: gfmat(nwhis,nwhis)
    real(8),allocatable:: frc(:),dfr(:),gfm(:)
    real(8):: ggg
    integer:: i,j
    allocate(frc(nwhis),dfr(nwhis),gfm(nwhis))
    do i=1,nwhis
      frc(i)=(frhis(i)+frhis(i+1))/2d0     ! bin-center frequency (Ha)
      dfr(i)= frhis(i+1)-frhis(i)          ! bin width (Ha)
    enddo
    do j=1,nwhis                            ! source bin (column)
       do i=1,nwhis                         ! target bin (row)
          gfm(i)= exp( -(frc(i)-frc(j))**2/(2d0*smearx0**2) ) * dfr(i)
       enddo
       ggg = sum(gfm(:))                    ! = sum_i K(i,j) dfr(i)
       do i=1,nwhis
          gfmat(i,j)= gfm(i)/ggg            ! = K(i,j) dfr(i) / sum_i' K(i',j) dfr(i')
       enddo
    enddo
    deallocate(frc,dfr,gfm)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!      do j=1,nwhis,20
! !      if(j>10) cycle
!       write(1019,*)
!       write(1019,*)
!       do i=1,nwhis
!         write(1019,'(2f19.8)') frhis(i),gfmat(i,j) !/(frhis(i+1)-frhis(i))
!       enddo
! !      write(*,*)'sssssss',i,sum(gfmat(:,j)*([(frhis(i+1)-frhis(i),i=1,nwhis)]))
!       write(*,*)'sssssss',j,sum(gfmat(:,j))
!     enddo
!     stop 'xxxxxxxxxxxaaa'
    
  end function gaussianfilterhis
end module m_dpsion


