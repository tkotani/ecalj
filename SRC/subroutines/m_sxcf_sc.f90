!! Calcualte full simga_ij(e_i)= <i|Re[Sigma](e_i)|j>  !checked 2023mar
!! ---------------------
!!   exchange=T: Calculate the exchange self-energy
!!   exchange=F: Calculate correlated part of the self-energy
!! \param zsec
!!   - S_ij= <i|Re[S](e_i)|j>
!!   - Note that S_ij itself is not Hermite becasue it includes e_i.
!!     i and j are band indexes
!! \remark
!!We now support only mode=3-----------------
!!  old version had 1,3,5scGW mode.
!!   diag+@EF      jobsw==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
!!   modeB (Not Available now)  jobsw==2 SE_nn'((e_n+e_n')/2)
!!   mode A        jobsw==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2 (Usually usued in QSGW).
!!   diagonly      jobsw==5 delta_nn' SE_nn(e_n) (not efficient memoryuse; but we don't use this mode so often).
!!
!!     zsec given in this routine is simply written as <i|Re[S](e_i)|j>.
!!     Be careful as for the difference between
!!     <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
!!     ---because e_i is included.
!!     The symmetrization (hermitian) procedure is inlucded in hqpe.sc.F
!!
!! Caution! npm=2 is not examined enough...
!!
!! Calculate the exchange part and the correlated part of self-energy.
!! T.Kotani started development after the analysis of F.Aryasetiawan's LMTO-ASA-GW.
!! We still use some of his ideas in this code.
!!
!! See paper
!! [1]T.Kotani, Quasiparticle Self-Consistent GW Method Based on the Augmented Plane-Wave
!!    and Muffin-Tin Orbital Method, J. Phys. Soc. Jpn., vol. 83, no. 9, p. 094711 [11 Pages], Sep. 2014.
!! [2]T.Kotani and M. van Schilfgaarde, Quasiparticle self-consistent GW method:
!!     A basis for the independent-particle approximation, Phys. Rev. B, vol. 76, no. 16, p. 165106[24pages], Oct. 2007.
!!=== Memo for Omega integral for SEc =====
!! The integral path is deformed along the imaginary-axis, but together with contribution of poles.
!!   See Fig.1 and around in Ref.[2].
!!1.Integration along the imaginary axis: Here is a memo originally by F.Aryasetiawan
!!     (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
!!  Gaussian integral along the imaginary axis.
!!  Transform: x = 1/(1+w')
!!    This leads to a denser mesh in w' around 0 for equal mesh x
!!    which is desirable since Wc and the lorentzian are peaked around w'=0
!!      wint = - (1/pi) < [x=0,1] Wc(iw') (w-e)x^2/{(w-e)^2 + w'^2} >
!!    The integrand is peaked around w'=0 or x=1 when w=e.
!!    To handel the problem, add and substract the singular part as follows:
!!     wint = - (1/pi) < [x=0,1] { Wc(iw') - Wc(0)exp(-a^2 w'^2) }
!!     * (w-e)/{(w-e)^2 +w'^2}x^2 > - (1/2) Wc(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|).
!!    The second term of the integral can be done analytically, which
!!    results in the last term a is some constant.
!!    When w = e, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w') and the integral becomes -Wc(0)/2
!!    This together with the contribution from the pole of G gives the so called static screened exchange -Wc(0).
!!2. Integration along real axis (contribution from the poles of G: SEc(pole))
!!    See Eq.(34),(55), and (58) and around in Ref.[2]. We now use Gaussian Smearing.
!
!!     q      = q-vector in SEc(q,t).
!!    itq     = ends states for SE
!!    ntq     = # of states t
!!    sxs_eq      = eigenvalues at q
!!     ef     = fermi level in Rydberg
!!   WVI, WVR: direct access files for W. along im axis (WVI) or along real axis (WVR)
!!   freq_r(nw_i:nw)   = frequencies along real axis. freq_r(0)=0d0
!!     wk     = weight for each k-point in the FBZ
!!    qbz     = k-points in the 1st BZ
!!     wx      = weights at gaussian points x between (0,1)
!!     ua_     = constant in exp(-ua^2 w'^2) s. wint.f
!!     expa    = exp(-ua^2 w'^2) s. wint.f
!!    irkip(k,R,nq) = gives index in the FBZ with k{IBZ, R=rotation
!!   nqibz   = number of k-points in the irreducible BZ
!!   nqbz    =                           full BZ
!!    nctot   = total no. of allowed core states
!!    nbloch  = total number of Bloch basis functions
!!    ndima   = total number of augmenation functions
!!    ngrp    = no. group elements (rotation matrices)
!!    niw     = no. frequencies along the imaginary axis
!!    nw_i:nw  = no. frequencies along the real axis. nw_i=0 or -nw.
!!    zsec(itp,itpp,iq)> = <psi(itp,q(:,iq)) |SEc| psi(iq,q(:,iq)>
!! \endverbatim
module m_sxcf_sc
  use m_readeigen, only: Readeval
  ! Step WA3/WB.3a: file I/O on __WVR.<kx> / __WVI.<kx> for reading goes through
  ! the m_wv_storage singleton.
  use m_wv_storage, only: wv_init_file, &
       wv_open_iq_for_read, wv_close_iq_for_read, wv_get_real, wv_get_imag
  use m_zmel, only: build_zmel, set_m2e_prod_basis, zmel, nbb !get_zmel_init =>
  use m_itq, only: ntq, nbandmx
  use m_struct_from_lmf,only: nspin,nband; use m_gw_product_basis,only: ndima; use m_core_state,only: nctot,ecore; use m_gw_user_config,only: niw
  use m_read_bzdata, only: qibz, qbz, wk=>wbz, nqibz, nqbz, wklm, lxklm, wqt=>wt
  use m_readVcoud, only: Readvcoud, ReleaseZcousq, vcoud, ngb, ngc
  use m_readfreq_r, only: freq_r, nw_i, nw, freqx, wx=>wwx, nblochpmx, mrecl, expa_, npm
  !use m_readhbe,only: nband
  use m_readgwinput, only: ua_, corehole, wcorehole
  use m_ftox
  use m_lgunit, only: stdo
  use m_sxcf_count,only: kxc, nstti, nstte, nstte2, nwxic, nwxc, icountini, icountend, irkip, ncount
  use m_nvfortran, only: findloc
  use m_hamindex, only: ngrp
  use m_blas, only: m_op_c, m_op_n, m_op_t, BACKEND_SIGMA
#if defined(__MP) && defined(__GPU)
  use m_blas, only: cmm_h16_d, sigma_fp16
#endif
use m_cmdopt_registry, only: c0_debug, c0_WVR2ptRaxis, c0_wcsmear
use m_GWinput, only: tg_wcsmear => wcsmear
use m_wfac, only: pole_weights
#if defined(__MP) && defined(__GPU)
  use m_blas, only: gemm => cmm_d
#elif defined(__MP)
  use m_blas, only: gemm => cmm_h
#elif defined(__GPU)
  use m_blas, only: gemm => zmm_d
#else
  use m_blas, only: gemm => zmm_h
#endif
  use m_kind, only: kp => kindgw
  !  use m_sxcf_main,only: zsecall
  use m_stopwatch
  use m_mem,only:writemem
  implicit none
  complex(kind=kp),allocatable,protected,target:: zsecall(:,:,:,:) !output
  public sxcf_scz_correlation, sxcf_scz_exchange, zsecall, reducez
  public sxcf_correlation_init, sxcf_correlation_step_kx, sxcf_correlation_finalize  ! Step WB.2a
  private
  real(8), parameter :: pi = 4d0*datan(1d0), fpi = 4d0*pi
  logical, parameter :: timemix = .true.
  complex(kind=kp), parameter :: CONE = (1_kp, 0_kp), CZERO = (0_kp, 0_kp)
  ! Loop iterators kx, irot, ip, isp are subroutine-local (declared inside
  ! each subroutine) so no time-dependent state leaks across subroutine
  ! boundaries.
  integer, allocatable :: ndiv(:), nstatei(:,:), nstatee(:,:)
  complex(kind=kp), allocatable :: wvr_upper(:,:), wvi_upper(:,:)
#ifdef __GPU
  attributes(device) :: wvr_upper, wvi_upper  ! GPU-only; no CPU copy
#endif
  ! Step WB.3b: working state lives at module scope (ecalj-style singleton).
  ! Lifetime: allocated in init / start of exchange, used by step_kx /
  ! kxloop body, deallocated in finalize / end of exchange.
  ! Stopwatches.
  type(stopwatch) :: sxs_zmel, sxs_xc, sxs_cr, sxs_ci, sxs_setwv, sxs_pole
  ! Shared workspace (used by both exchange and correlation flows).
  real(8), allocatable :: sxs_ekc(:), sxs_eq(:)
  integer :: sxs_ntqxx
  real(8) :: sxs_wkkr
  ! Correlation-only workspace.
  real(8), allocatable :: sxs_omega(:)
  logical :: sxs_keepwv
  integer :: sxs_nt0p, sxs_nt0m
  ! Omega-mesh split bounds (computed once in _init, reused per-kx in _step_kx).
  integer :: sxs_wi_ini = 0, sxs_wi_fin = 0, sxs_wi_num = 0
  integer :: sxs_wr_ini = 0, sxs_wr_fin = 0, sxs_wr_num = 0
  real(8), parameter:: rmax=2d0
contains
  subroutine reducez(nspinmx)
#ifdef __MP
    use m_mpi,only: MPI__reduceSum => MPI__reduceSum_c
#else
    use m_mpi,only: MPI__reduceSum
#endif
    integer::nspinmx
    call MPI__reduceSum(root=0, data=zsecall, sizex=ntq*ntq*nqibz*nspinmx )
  end subroutine reducez
  subroutine sxcf_scz_exchange(ef, esmr, ixc, nspinmx) !ixc is dummy
    implicit none
    integer :: icount, ns1, ns2, kr, izz
    integer :: kx, irot, ip, isp  ! Step 2.1: localized from module scope
    logical, parameter :: debug=.false.
    integer, intent(in) :: nspinmx, ixc
    real(8), intent(in) :: ef, esmr
    real(8) :: q(3), qibz_k(3), qbz_kr(3), qk(3)
    character(64):: charli
    character(8):: charext
    allocate(sxs_ekc(nctot+nband), sxs_eq(nband))
    if(nw_i/=0) call rx('Current version we assume nw_i=0. Time-reversal symmetry')
    LoopScheduleCheck: block
      izz=0
      kxloopX:                do kx  =1,nqibz
        irotloopX:            do irot=1,ngrp
          iploopexternalX:    do ip=1,nqibz
            isploopexternalX: do isp=1,nspinmx
              kr = irkip(isp,kx,irot,ip)
              if(kr==0) cycle
              NMBATCHloopX:   do icount = icountini(isp,ip,irot,kx),icountend(isp,ip,irot,kx) !batch of middle states.
                izz=izz+1
                write(stdo,ftox)'= kxloop Schedule ',izz, ' iqibz irot ip isp icount=',kx,irot,ip,isp,icount
              enddo NMBATCHloopX
            enddo isploopexternalX
          enddo iploopexternalX
        enddo irotloopX
      enddo kxloopX
    end block LoopScheduleCheck
    izz=0
    call stopwatch_init(sxs_zmel, 'zmel')
    call stopwatch_init(sxs_xc, 'ex')
    if (allocated(zsecall)) then
       !$acc exit data delete(zsecall)
       deallocate(zsecall)
    endif
    allocate(zsecall(ntq,ntq,nqibz,nspinmx))
    !$acc enter data create(zsecall)
    !$acc kernels
    zsecall(1:ntq,1:ntq,1:nqibz,1:nspinmx) = CZERO
    !$acc end kernels
    kxloop:                do kx=1, nqibz      ! kx is irreducible !kx is main axis where we calculate W(kx).
      qibz_k = qibz(:,kx)
      call Readvcoud(qibz_k, kx, NoVcou=.false.)   !Readin ngc,ngb,vcoud ! Coulomb matrix
      call set_m2e_prod_basis(npr=ngb)     !Set M to E basis transformation matrix
      call ReleaseZcousq()                         !Release zcousq used in set_m2e_prod_basis
      irotloop:            do irot=1, ngrp     ! (kx,irot) determines qbz(:,kr), which is in FBZ. W(kx) is rotated to be W(g(kx))
        iploopexternal:    do ip=1, nqibz      !external index for q of \Sigma(q,isp)
          isploopexternal: do isp=1, nspinmx   !external index
            kr = irkip(isp,kx,irot,ip)
            if(kr==0) cycle
            q = qibz(:,ip)
            qbz_kr = qbz (:,kr)   !rotated qbz vector.
            qk = q - qbz_kr       !<M(qbz_kr) phi(q-qbz_kr)|phi(q)>
            sxs_eq = readeval(q,isp)  !readin eigenvalue
            sxs_ekc(1:nctot+nband) = [ecore(1:nctot,isp),readeval(qk, isp)]
            sxs_ntqxx = nbandmx(ip,isp) ! sxs_ntqxx is number of bands for <i|sigma|j>.
            sxs_wkkr = wk(kr)
            NMBATCHloop:   do icount = icountini(isp,ip,irot,kx), icountend(isp,ip,irot,kx) !batch of middle states.
              ns1 = nstti(icount)  !Range of middle states is [ns1:ns2] for given icount
              ns2 = nstte(icount)  !
              call stopwatch_start(sxs_zmel)
              izz=izz+1
              call writemem('=== KXloop '//trim(charext(izz))//' iqiqz irot ip isp icount= '//&
                   trim(charli([kx,irot,ip,isp,icount],5)))
              call build_zmel(q,qibz_k,irot,qbz_kr,ns1,ns2,isp,1,sxs_ntqxx,isp,nctot,ncc=0,zmelconjg=.false., &
                                      is_m_basis=.false., mpi_mode=.false.)
              call writemem('    endof build_zmel')
              call stopwatch_pause(sxs_zmel)
              call stopwatch_start(sxs_xc)
              associate( zsec=>zsecall(:,:,ip,isp) )
                get_exchange_block:block !subroutine get_exchange(ef, esmr, ns1, ns2, zsec)
                  real(8) :: wfacx, wtff(ns1:ns2) ! external function
                  real(8), allocatable :: vcoud_buf(:)
                  complex(kind=kp), allocatable :: vzmel(:,:,:)
                  integer :: it, itp, itpp, ierr,is1
#ifdef __GPU
                  attributes(device) :: vzmel, vcoud_buf
#endif
                  if(ns1 > ns2) goto 1110 !instead of return. Use guard clause coding.
                  do is1=ns1,ns2
                    if(is1<=nctot) then; wtff(is1) = 1d0 !these are for nvfortran24.1
                    else;                wtff(is1) = wfacx(-1d99, ef, sxs_ekc(is1+nctot), esmr)
                    endif
                  enddo
                  if(corehole) wtff(ns1:nctot) = wtff(ns1:nctot) * wcorehole(ns1:nctot,isp)
                  allocate(vcoud_buf(ngb))
                  !$acc data copyin(vcoud, wklm(1), wk(1), wtff) present(zmel)
                  !$acc kernels
                  vcoud_buf(1:ngb) = vcoud(1:ngb)
                  !$acc end kernels
                  if(kx == 1) vcoud_buf(1) = wklm(1)*fpi*sqrt(fpi)/wk(1) ! voud_buf(1) is effective v(q=0) in the Gamma cell.
                  allocate(vzmel(1:nbb,ns1:ns2,1:sxs_ntqxx))
                  !$acc kernels loop independent collapse(2)
                  do itpp = 1, sxs_ntqxx
                    do it = ns1, ns2
                      vzmel(1:ngb,it,itpp) = cmplx(wtff(it)*vcoud_buf(1:ngb)*zmel(1:ngb,it,itpp),kind=kp)
                    enddo
                  enddo
                  !$acc end kernels
                  !$acc host_data use_device(zmel, zsec)
                  ierr = gemm(zmel, vzmel, zsec, sxs_ntqxx, sxs_ntqxx, (ns2-ns1+1)*nbb, opA = m_op_c, &
                       alpha = cmplx(-sxs_wkkr,0_kp,kind=kp), beta = CONE, ldC = ntq, splitk = ns2-ns1+1) !one k range per middle state (nbb each)
                  !$acc end host_data
                  !$acc end data
                  deallocate(vzmel, vcoud_buf)!, wtff)
1110              continue
                end block get_exchange_block !end subroutine get_exchange
              endassociate
              call writemem('    endof ExchangeSelfEnergy')
              call stopwatch_pause(sxs_xc)
              write(stdo,ftox) '    End of icount:', icount ,' of', ncount, &
                   'zmel:', ftof(stopwatch_lap_time(sxs_zmel),4), '(sec)', &
                   'exch:', ftof(stopwatch_lap_time(sxs_xc),4),   '(sec)'
              call flush(stdo)
            enddo NMBATCHloop
          enddo isploopexternal
        enddo iploopexternal
      enddo irotloop
    enddo kxloop
    !$acc exit data copyout (zsecall)
    deallocate(sxs_ekc, sxs_eq)
    call stopwatch_show(sxs_zmel)
    call stopwatch_show(sxs_xc)
  endsubroutine sxcf_scz_exchange

  ! ============================================================================
  ! Step WB.2a: sxcf_scz_correlation is now a thin wrapper around the
  ! _init / _step_kx / _finalize triple. The original behavior is preserved
  ! exactly: caller-side iteration over kx is identical to the old internal
  ! kxloop, and state lifetime matches the old subroutine-local lifetime.
  !
  ! The split exists so that streaming flows (Phase 1-C) can interleave
  ! per-iq W(0,kx) production with sxcf consumption: caller produces W at kx,
  ! calls step_kx(kx), and discards W before the next kx — without ever
  ! buffering the full per-iq W in __WVR/__WVI files or in 4D arrays.
  ! ============================================================================
  subroutine sxcf_scz_correlation(ef, esmr, ixc, nspinmx)
    integer, intent(in) :: nspinmx, ixc
    real(8), intent(in) :: ef, esmr
    integer :: kx
    call wv_init_file(mreclx=mrecl, nw_i=nw_i)
    call sxcf_correlation_init(ef, esmr, nspinmx)
    do kx = 1, nqibz
      call sxcf_correlation_step_kx(kx, ef, esmr, nspinmx)
    enddo
    call sxcf_correlation_finalize()
  end subroutine sxcf_scz_correlation

  ! Pre-kxloop setup: allocate module-level workspace, initialize stopwatches,
  ! partition omega-mesh across MPI w-ranks, allocate + zero zsecall.
  subroutine sxcf_correlation_init(ef, esmr, nspinmx)
    use m_keyvalue, only: getkeyvalue
    use m_GWinput, only: gwinput_init, gwinput_loaded, tg_KeepWV => KeepWV
    use m_mpi, only: mpi__size_b => mpi__size_b_sxc, mpi__rank_b => mpi__rank_b_sxc, ipr
    use m_blas, only: int_split
    use m_gpu, only: use_gpu
    real(8), intent(in) :: ef, esmr
    integer, intent(in) :: nspinmx
    integer :: kx, irot, ip, isp, kr, izz, icount
    if (nw_i /= 0) call rx('Current version we assume nw_i=0. Time-reversal symmetry')
    allocate(sxs_ekc(nctot+nband), sxs_eq(nband), sxs_omega(ntq))
    !$acc enter data create(sxs_ekc, sxs_omega) copyin(freqx, wx, expa_)
    call gwinput_init()
    if (gwinput_loaded) then
       sxs_keepwv = tg_KeepWV
    else
       call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
    endif
    LoopScheduleCheck: block
      izz = 0
      kxloopX:                do kx   = 1, nqibz
        irotloopX:            do irot = 1, ngrp
          iploopexternalX:    do ip   = 1, nqibz
            isploopexternalX: do isp  = 1, nspinmx
              kr = irkip(isp,kx,irot,ip)
              if (kr == 0) cycle
              NMBATCHloopX: do icount = icountini(isp,ip,irot,kx), icountend(isp,ip,irot,kx)
                izz = izz + 1
                write(stdo,ftox) '= KXloop Scheduling ', izz, ' iqiqz irot ip isp icount=', kx, irot, ip, isp, icount
              enddo NMBATCHloopX
            enddo isploopexternalX
          enddo iploopexternalX
        enddo irotloopX
      enddo kxloopX
    end block LoopScheduleCheck
    call int_split(    niw+1, mpi__size_b, mpi__rank_b, sxs_wi_ini, sxs_wi_fin, sxs_wi_num, start_index=0)
    call int_split(nw-nw_i+1, mpi__size_b, mpi__rank_b, sxs_wr_ini, sxs_wr_fin, sxs_wr_num, start_index=nw_i)
    if (ipr) write(stdo,ftox) 'Imag sxs_omega mesh split:', sxs_wi_ini, sxs_wi_fin, sxs_wi_num, &
         'Real sxs_omega mesh split:', sxs_wr_ini, sxs_wr_fin, sxs_wr_num
    if (ipr) write(stdo,ftox) '# of tasks:', izz
    call flush(stdo)
    call stopwatch_init(sxs_zmel,  'zmel')
    call stopwatch_init(sxs_xc,    'ec')
    call stopwatch_init(sxs_cr,    'ec realaxis (host)', hostonly=.true.)   ! inside the async batch: host time only;
    call stopwatch_init(sxs_ci,    'ec imagaxis (host)', hostonly=.true.)   ! 'ec' (sxs_xc) is the batch with the GPU
    call stopwatch_init(sxs_setwv, 'read wv')
    call stopwatch_init(sxs_pole,  'realaxis pole weights (CPU)', hostonly=.true.)
    if (allocated(zsecall)) then
       !$acc exit data delete(zsecall)
       deallocate(zsecall)
    endif
    allocate(zsecall(ntq,ntq,nqibz,nspinmx))
    !$acc enter data create(zsecall)
    !$acc kernels
    zsecall(1:ntq,1:ntq,1:nqibz,1:nspinmx) = CZERO
    !$acc end kernels
  end subroutine sxcf_correlation_init

  ! One iteration of the kxloop: read/load W(kx), then accumulate the
  ! correlation contribution into zsecall(:,:,ip,isp) for all (irot, ip, isp).
  subroutine sxcf_correlation_step_kx(kx, ef, esmr, nspinmx)
    use m_mpi, only: comm_b => comm_b_sxc, ipr
    use m_gpu, only: use_gpu
    use m_wfac, only: sig_window
    integer, intent(in) :: kx, nspinmx
    real(8), intent(in) :: ef, esmr
    integer :: icount, ns1, ns2, kr, nwxi, ns2r, nwx, izz, n_nttp, tri_idx
    integer :: irot, ip, isp
    real(8) :: q(3), qibz_k(3), qbz_kr(3), qk(3)
    logical :: debug
    integer, allocatable :: idx_i(:), idx_j(:)
    character(64) :: charli
    character(8)  :: charext
    debug = c0_debug
    qibz_k = qibz(:,kx)
    call Readvcoud(qibz_k, kx, NoVcou=.false.)   !Readin ngc,ngb,vcoud ! Coulomb matrix
    call set_m2e_prod_basis(npr=ngb)             !Set M to E basis transformation matrix
    call ReleaseZcousq()                         !Release zcousq used in set_m2e_prod_basis
    SetWVblock: block !subroutine setwv()
      integer :: iqini, iqend, iw, i, j
      character(10) :: i2char
      complex(kind=kp), allocatable :: wv(:,:)
      real(8), parameter :: gb = 1000*1000*1000
      if (any(kx == kxc(:))) then
        call wv_open_iq_for_read(kx, want_real=.true., want_imag=.true.)
        if (sxs_keepwv) then
          if (ipr) write(stdo,ftox) 'save WVI and WVR on GPU device memory (slice-by-slice from SHM, no CPU copy)'
          call flush(stdo)
          call stopwatch_reset(sxs_setwv)
          call stopwatch_start(sxs_setwv)
          allocate(idx_i(ngb*(ngb+1)/2), idx_j(ngb*(ngb+1)/2))
          tri_idx = 1
          do j = 1, ngb
            do i = 1, j
              idx_i(tri_idx) = i
              idx_j(tri_idx) = j
              tri_idx = tri_idx + 1
            enddo
          enddo
          !$acc enter data copyin(idx_i, idx_j)
          ! Allocate wvi/wvr_upper on GPU device memory only (attributes(device)),
          ! then fill slice-by-slice via small CPU buffer wv(ngb,ngb).
          allocate(wvi_upper(ngb*(ngb+1)/2, sxs_wi_ini:sxs_wi_fin))
          allocate(wvr_upper(ngb*(ngb+1)/2, sxs_wr_ini:sxs_wr_fin))
          allocate(wv(nblochpmx,nblochpmx))
          !$acc enter data create(wv)
          do iw = sxs_wi_ini, sxs_wi_fin
            if (iw == 0) then
              call wv_get_real(iw, wv)
            else
              call wv_get_imag(iw, wv)
            endif
            !$acc update device(wv)
            !$acc parallel loop present(wvi_upper, idx_i, idx_j, wv)
            do tri_idx = 1, ngb*(ngb+1)/2
              wvi_upper(tri_idx,iw) = wv(idx_i(tri_idx), idx_j(tri_idx))
            enddo
            !$acc end parallel loop
          enddo
          do iw = sxs_wr_ini, sxs_wr_fin
            call wv_get_real(iw, wv)
            !$acc update device(wv)
            !$acc parallel loop present(wvr_upper, idx_i, idx_j, wv)
            do tri_idx = 1, ngb*(ngb+1)/2
              wvr_upper(tri_idx,iw) = (wv(idx_i(tri_idx),idx_j(tri_idx)) + &
                                        conjg(wv(idx_j(tri_idx),idx_i(tri_idx)))) * 0.5_kp
            enddo
            !$acc end parallel loop
          enddo
          !$acc exit data delete(wv)
          deallocate(wv)
          call stopwatch_pause(sxs_setwv)
          if (ipr) write(stdo, '(X,A,2F8.3)') 'WVI/WVR GPU sizes (GB)', &
            dble(size(wvi_upper))*kp*2/gb, dble(size(wvr_upper))*kp*2/gb
          call stopwatch_show(sxs_setwv)
        endif
      endif
      !end subroutine setwv
    end block SetWVblock
    izz = 0
    ! --skipq0Sc (diagnostic): leave out the Gamma-cell W (kx=1, whose head comes from the offset-Gamma
    ! chi0) from Sigma_c.  Everything else (WV open/close, exchange) is untouched.
    irotloop:            do irot = 1, ngrp    ! (kx,irot) determines qbz(:,kr), which is in FBZ. W(kx) is rotated to be W(g(kx))
      iploopexternal:    do ip   = 1, nqibz   !external index for q of \Sigma(q,isp)
        isploopexternal: do isp  = 1, nspinmx !external index
          kr = irkip(isp,kx,irot,ip)
          if (kr == 0) cycle
          q = qibz(:,ip)
          qbz_kr = qbz(:,kr)   !rotated qbz vector.
          qk = q - qbz_kr      !<M(qbz_kr) phi(q-qbz_kr)|phi(q)>
          sxs_eq = readeval(q,isp)  !readin eigenvalue
          sxs_wkkr = wk(kr)
          sxs_ekc(1:nctot+nband) = [ecore(1:nctot,isp), readeval(qk, isp)]
          sxs_ntqxx = nbandmx(ip,isp) ! sxs_ntqxx is number of bands for <i|sigma|j>.
          sxs_omega(1:ntq) = sxs_eq(1:ntq)
          !$acc update device(sxs_ekc, sxs_omega)
          sxs_nt0p = count(sxs_ekc < ef + sig_window(esmr))   ! states that can be partially occupied
          sxs_nt0m = count(sxs_ekc < ef - sig_window(esmr))   ! (kernel tail, 15 kBT of t_sigmaw)
          NMBATCHloop: do icount = icountini(isp,ip,irot,kx), icountend(isp,ip,irot,kx) !batch of middle states.
            ns1  = nstti(icount)   !Range of middle states is [ns1:ns2] for given icount
            ns2  = nstte(icount)
            nwxi = nwxic(icount)   !minimum sxs_omega for W
            nwx  = nwxc(icount)    !max sxs_omega for W
            ns2r = nstte2(icount)  !Range of middle states [ns1:ns2r] for CorrelationSelfEnergyRealAxis
            izz  = izz + 1
            call writemem('=== KXloop '//trim(charext(izz))//' iqiqz irot ip isp icount= '//&
                 trim(charli([kx,irot,ip,isp,icount],5)))
            call stopwatch_start(sxs_zmel)
            call build_zmel(q,qibz_k,irot,qbz_kr,ns1,ns2,isp,1,sxs_ntqxx,isp,nctot,ncc=0,zmelconjg=.false., &
                                 is_m_basis=.true., mpi_mode=.not.use_gpu, comm=comm_b)
            call writemem('    endof build_zmel')
            call stopwatch_pause(sxs_zmel)
            call stopwatch_reset(sxs_setwv)
            call stopwatch_start(sxs_xc)
            ! call get_correlation(ef, esmr, ns1, ns2, ns2r, nwxi, nwx, zsecall(1,1,ip,isp))
            associate( zsec=>zsecall(:,:,ip,isp) )
              get_correlation_block :block
                ! On the GPU the whole batch runs on OpenACC queue 1 (sigma_stream_begin): the kernels are async(1) and
                ! the products of m_blas go to the same stream, so they stay in order without a host wait per kernel.
                ! While the imaginary-axis products run, the host computes the pole weights of the real axis; no device
                ! allocation or free in that window (it would wait for the device).
                real(8), parameter :: wfaccut=1d-8
                complex(kind=kp), parameter :: img=(0_kp,1_kp)
                complex(kind=kp) :: beta
                complex(kind=kp), allocatable :: czmelwc(:,:,:), wzmel(:,:,:), wz_iw(:,:), czwc_iw(:,:)
                integer :: it, itp, iw, ierr, i, j, nttp_max, nttp(0:nw), igb, ntw
#if defined(__MP) && defined(__GPU)
                ! FP16 route of the Sigma_c products (realhgemm): the weighting kernels write B in FP16 themselves, scaled
                ! by a power of 2 from the bound max|w| max|zmel| (no scan of B, no conversion kernel).  cmm_h16_d.
                real(2), allocatable :: bh(:)
                attributes(device) :: bh
                real(8) :: zmx, wmx
                real(8), allocatable :: wmaxiw(:)
                real(4) :: fbw
                logical :: h16
                complex(4) :: v16
                integer :: j16
#endif
                complex(kind=kp), allocatable :: wv(:,:), wc(:,:)
                real(8), allocatable :: wgtim(:,:,:), wgtiw(:,:)
                integer, allocatable :: itw(:,:), itpw(:,:)
                real(kind=kp) :: zsec_img
#ifdef __GPU
                attributes(device) :: czmelwc, wc, wzmel, wz_iw, czwc_iw
#endif
                if (ns1 > ns2) goto 1114 !instead of return
                allocate(wv(nblochpmx,nblochpmx))
                allocate(wc(ngb,ngb))
                allocate(czmelwc, mold = zmel)
                allocate(wzmel(1:ngb,ns1:ns2,1:sxs_ntqxx))
                allocate(wgtim(0:npm*niw,ns1:ns2,sxs_ntqxx))
                !$acc enter data create(wgtim)
                nttp = 0
                nttp_max = 0
                call sigma_stream_begin()
#if defined(__MP) && defined(__GPU)
                allocate(bh(2*ngb*(ns2-ns1+1)*sxs_ntqxx), wmaxiw(0:npm*niw))
                zmx = 0d0
                !$acc parallel loop collapse(3) present(zmel) reduction(max:zmx)
                do itp = 1, sxs_ntqxx
                  do it = ns1, ns2
                    do igb = 1, ngb
                      zmx = max(zmx, dble(abs(real(zmel(igb,it,itp)))), dble(abs(aimag(zmel(igb,it,itp)))))
                    enddo
                  enddo
                enddo
#endif
                call stopwatch_start(sxs_ci)
                CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)! Integration along ImAxis for zwz(sxs_omega)
                  use m_readfreq_r, only: wt=>wwx, x=>freqx
                  use m_wfac, only: fd_cdf
                  real(8) :: we, aw, aw2, u, xk, estep, cons, w1, s0, omg, ek, wkkr, uaa, esm
                  integer :: jq, nq, ntqxx, niwx, npmx, nct, ks1, ks2
                  ! Imaginary-axis weights of one level ekc(it) for Sigma_c(omega(itp)), Eq. 57 of PRB 76, 165106:
                  ! numeric part on the niw mesh with W_c(i w') - W_c(0) exp(-(ua w')^2) (smooth), analytic
                  ! part for the Gaussian fit W_c(0) exp(-(ua w')^2): -sign(we)/2 exp(aw^2) erfc(aw), aw=ua|we|,
                  ! we=(omega-e)/2 (Hartree).  This is the sharp-level formula (it was used for the core levels;
                  ! the valence levels had a Gaussian regularization 'sig' that smeared the sign(we) step
                  ! and was the Gaussian level smearing of the esmr era).
                  ! The level smearing is now the Fermi-Dirac kernel of width esmr (=kBT, m_wfac), the same
                  ! kernel as the pole-term weights and wcsmear: the weights of a valence level within
                  ! 30 kBT of omega are averaged over the kernel, e = ekc + x, x = kBT ln(u/(1-u)) with
                  ! midpoint nodes u in (0,1) (the FD cumulative is the measure).  The step -sign(we)/2 is
                  ! discontinuous and is averaged analytically, -(2 Phi_FD(omega-ekc) - 1)/2; the smooth rest
                  ! -sign(we)/2 (exp(aw^2) erfc(aw) - 1) and the numeric part go through the nodes.  So
                  ! imaginary-axis + pole term = Sigma_c of the FD-smeared level (with wcsmear), and the two
                  ! half-residue steps cancel for a level near omega (NiO 2^3 O 2s pair, 2026-09-20).
                  ! Mechanism: for one sharp level, I(we) (this block) has the step -sign(we) W_c(0)/2 at we=0
                  ! and the pole term P has the step 0 -> W_c(0) at the window edge e = omega; the two cancel.
                  ! Smearing the level = averaging BOTH over the kernel; a kernel change here must be the same
                  ! as in wfacx2/pole_weights (contract in m_wfac, wfacx.f90).  Formulae and the NiO numbers:
                  ! https://ecalj.github.io/ecaljdoc/manual/kBT#_3-6-sigma-c-の-contour-分解と準位-smearing-の整合
                  ! One (it,itp) per GPU thread; on the CPU the same loops run in order.
                  integer, parameter :: nqfd = 40
                  ntqxx = sxs_ntqxx; niwx = niw; npmx = npm; nct = nctot; ks1 = ns1; ks2 = ns2
                  wkkr = sxs_wkkr; uaa = ua_; esm = esmr
                  !$acc parallel loop gang vector collapse(2) async(1) present(wgtim, x, wt, expa_, sxs_omega, sxs_ekc) &
                  !$acc   private(we, aw, aw2, u, xk, estep, cons, w1, s0, omg, ek, jq, nq, iw)
                  itpdo: do itp = 1, ntqxx
                    itpo: do it = ks1, ks2
                      omg = sxs_omega(itp)
                      ek  = sxs_ekc(it)
                      nq = 1
                      if (it>nct .and. esm>0d0 .and. abs(omg-ek) < 30d0*esm) nq = nqfd
                      !$acc loop seq
                      do iw = 0, npmx*niwx
                        wgtim(iw,it,itp) = 0d0
                      enddo
                      !$acc loop seq
                      fdnodes: do jq = 1, nq
                        xk = 0d0
                        if (nq > 1) then
                          u  = (jq - .5d0)/nq
                          xk = esm*log(u/(1d0-u))
                        endif
                        we = .5d0*(omg - ek - xk) !we in hartree unit (atomic unit)
                        aw = abs(uaa*we)
                        aw2 = aw*aw
                        s0 = 0d0
                        !$acc loop seq
                        do iw = 1, niwx
                          cons = 1d0/(we**2*x(iw)**2 + (1d0-x(iw))**2)          ! = 1/(x^2 (w'^2+we^2)), w' = 1/x - 1
                          w1 = we*cons*wt(iw)*(-1d0/pi)
                          wgtim(iw,it,itp) = wgtim(iw,it,itp) + w1/nq
                          s0 = s0 + w1*expa_(iw)
                          if (npmx==2) wgtim(niwx+iw,it,itp) = wgtim(niwx+iw,it,itp) + cons*(1d0/x(iw)-1d0)*wt(iw)/pi/nq !Asymmetric contribution need check
                        enddo
                        estep = dsign(1d0,we)*dexp(aw2)*erfc(aw)     ! sign(we) exp(aw^2) erfc(aw), -> sign(we) at we=0
                        if (nq > 1) estep = estep - dsign(1d0,we)    ! the step itself is added analytically below
                        if (dabs(we) < rmax/uaa) wgtim(0,it,itp) = wgtim(0,it,itp) + (-s0 - 0.5d0*estep)/nq
                      enddo fdnodes
                      if (nq > 1) wgtim(0,it,itp) = wgtim(0,it,itp) - 0.5d0*(2d0*fd_cdf(omg-ek, esm) - 1d0)
                      !$acc loop seq
                      do iw = 0, npmx*niwx
                        wgtim(iw,it,itp) = wkkr*wgtim(iw,it,itp) !! Integration weight wgtim along im axis for zwz(0:niw*npm)
                      enddo
                    enddo itpo
                  enddo itpdo
#if defined(__MP) && defined(__GPU)
                  !$acc wait(1)
                  !$acc parallel loop gang present(wgtim) copyout(wmaxiw(0:npmx*niwx)) private(wmx)
                  do iw = 0, npmx*niwx
                    wmx = 0d0
                    !$acc loop vector collapse(2) reduction(max:wmx)
                    do itp = 1, ntqxx
                      do it = ks1, ks2
                        wmx = max(wmx, abs(wgtim(iw,it,itp)))
                      enddo
                    enddo
                    wmaxiw(iw) = wmx
                  enddo
#endif
                  if (debug) call writemem('    Goto iwimag')
                  if (debug) write(stdo,ftox) 'mmmmSc size of mm in imagaxis', (ns2-ns1+1)*sxs_ntqxx, ngb, ngb
                  iwimag: do iw = sxs_wi_ini, sxs_wi_fin ! iwimag:do iw = 0, niw !niw is ~10. ixx=0 is for sxs_omega=0 nw_i=0 (Time reversal) or nw_i =-nw
                    if (iw < 0 .or. iw > niw) cycle
                    if (sxs_keepwv) then
                      !$acc kernels loop independent present(wvi_upper, idx_i, idx_j) async(1)
                      do tri_idx = 1, ngb*(ngb+1)/2
                        i = idx_i(tri_idx)
                        j = idx_j(tri_idx)
                        wc(i,j) = wvi_upper(tri_idx,iw)
                        wc(j,i) = conjg(wvi_upper(tri_idx,iw))
                      enddo
                      !$acc end kernels
                    else
                      !$acc wait(1)
                      call stopwatch_start(sxs_setwv)
                      if (iw == 0) call wv_get_real(iw, wv)
                      if (iw > 0)  call wv_get_imag(iw, wv)
                      wc(1:ngb,1:ngb) = wv(1:ngb,1:ngb)  !copy to GPU
                      call stopwatch_pause(sxs_setwv)
                    endif
                    beta = CONE
                    if (iw == sxs_wi_ini) beta = CZERO
#if defined(__MP) && defined(__GPU)
                    h16 = sigma_fp16(ngb, (ns2-ns1+1)*sxs_ntqxx, ngb)
                    if (h16) then
                      fbw = pow2scale(wmaxiw(iw)*zmx)
                      !$acc parallel loop collapse(3) present(zmel, wgtim) private(v16, j16) async(1)
                      do itp = 1, sxs_ntqxx
                        do it = ns1, ns2
                          do igb = 1, ngb
                            v16 = cmplx(real(wgtim(iw,it,itp),4)*fbw*zmel(igb,it,itp), kind=4)
                            j16 = 2*(igb + ngb*((it-ns1) + (ns2-ns1+1)*(itp-1))) - 1
                            bh(j16)   = real(real(v16), kind=2)
                            bh(j16+1) = real(aimag(v16), kind=2)
                          enddo
                        enddo
                      enddo
                      ierr = cmm_h16_d(wc, bh, fbw, czmelwc, ngb, (ns2-ns1+1)*sxs_ntqxx, ngb, opa=m_op_C, beta=beta, &
                                       key = 1000 + iw)
                      cycle
                    endif
#endif
                    !$acc parallel loop collapse(3) present(zmel, wgtim) async(1)
                    do itp = 1, sxs_ntqxx
                      do it = ns1, ns2
                        do igb = 1, ngb
                          wzmel(igb,it,itp) = cmplx(wgtim(iw,it,itp)*zmel(igb,it,itp), kind=kp)
                        enddo
                      enddo
                    enddo
                    !the most time-consuming part in the correlation part
                    ierr = gemm(wc, wzmel, czmelwc, ngb, (ns2-ns1+1)*sxs_ntqxx, ngb, beta = beta, opA = m_op_C, &
                                key = 1000 + iw, policy = BACKEND_SIGMA)   ! W(i omega) fixed for this kx (m_zmel resets keys)
                  enddo iwimag
                EndBlock CorrelationSelfEnergyImagAxis
                if (debug) call writemem('    endof CorrelationSelfEnergyImagAxis')
                call stopwatch_pause(sxs_ci)

                call stopwatch_start(sxs_cr)
                CorrelationSelfEnergyRealAxis: Block !Real Axis integral. Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
                  use m_wfac, only: wfacx2
                  integer :: itini, itend, ittp, i, iw1, iw2, ipass
                  real(8) :: omg, wfac, wts(0:nw), esmr_it
                  logical :: smear
                  smear = tg_wcsmear .or. c0_wcsmear
                  ! On the GPU the imaginary-axis products above are still running on queue 1 while the host
                  ! computes these weights.
                  call stopwatch_start(sxs_pole)
                  PoleWeights: do ipass = 1, 2   ! pass 1 counts the pairs per mesh point, pass 2 stores them
                    nttp = 0
                    do itp = 1, sxs_ntqxx
                      omg   = sxs_omega(itp)
                      itini = merge(max(ns1,sxs_nt0m+1),  ns1, mask= omg>=ef)
                      itend = merge(ns2r,  min(sxs_nt0p,ns2r), mask= omg>=ef)
                      do it = itini, itend     ! sxs_nt0p corresponds to efp
                        esmr_it = merge(0d0, esmr, mask=it<=nctot)   ! core levels are sharp
                        wfac = wfacx2(omg, ef, sxs_ekc(it), esmr_it)  ! weight of the level inside [ef, omg]
                        if (wfac < wfaccut) cycle
                        call pole_weights(omg, ef, sxs_ekc(it), esmr_it, smear, nw, freq_r(0:nw), wfac, &
                                          sxs_wkkr*dsign(1d0, omg-ef), iw1, iw2, wts)
                        do i = iw1, iw2
                          if (wts(i) == 0d0) cycle
                          nttp(i) = nttp(i) + 1
                          if (ipass == 2) then
                            itw(nttp(i),i)   = it
                            itpw(nttp(i),i)  = itp
                            wgtiw(nttp(i),i) = wts(i)
                          endif
                        enddo
                      enddo
                    enddo
                    if (ipass == 1) then
                      nttp_max = maxval(nttp)
                      if (nttp_max <= 0) then
                        call stopwatch_pause(sxs_pole)
                        goto 1113
                      endif
                      allocate (itw(nttp_max,0:nw),  source = 0)
                      allocate (itpw(nttp_max,0:nw), source = 0)
                      allocate (wgtiw(nttp_max,0:nw),source = 0d0)
                    endif
                  enddo PoleWeights
                  call stopwatch_pause(sxs_pole)
                  n_nttp = count(nttp(sxs_wr_ini:sxs_wr_fin) > 0)
                  allocate(wz_iw(ngb,nttp_max), czwc_iw(ngb,nttp_max))
                  !$acc enter data copyin(wgtiw, nttp, itw, itpw)
                  iwreal: do iw = sxs_wr_ini, sxs_wr_fin
                    if (iw < nwxi .or. iw > nwx) cycle
                    if (nttp(iw) < 1) cycle
                    if (sxs_keepwv) then
                      !$acc kernels loop independent present(wvr_upper, idx_i, idx_j) async(1)
                      do tri_idx = 1, ngb*(ngb+1)/2
                        i = idx_i(tri_idx)
                        j = idx_j(tri_idx)
                        wc(i,j) = wvr_upper(tri_idx,iw)
                        wc(j,i) = conjg(wvr_upper(tri_idx,iw))
                      enddo
                      !$acc end kernels
                    else
                      !$acc wait(1)
                      call stopwatch_start(sxs_setwv)
                      call wv_get_real(iw, wv)
                      wc(1:ngb,1:ngb) = wv(1:ngb,1:ngb)  !copy to GPU
                      !$acc kernels
                      wc(:,:) = (wc(:,:) + transpose(conjg(wc(:,:))))*0.5_kp
                      !$acc end kernels
                      call stopwatch_pause(sxs_setwv)
                    endif
                    ntw = nttp(iw)
#if defined(__MP) && defined(__GPU)
                    h16 = sigma_fp16(ngb, ntw, ngb)
                    if (h16) then
                      fbw = pow2scale(maxval(abs(wgtiw(1:ntw,iw)))*zmx)
                      !$acc parallel loop collapse(2) present(zmel, wgtiw, itw, itpw) private(v16, j16) async(1)
                      do ittp = 1, ntw
                        do igb = 1, ngb
                          v16 = cmplx(real(wgtiw(ittp,iw),4)*fbw*zmel(igb,itw(ittp,iw),itpw(ittp,iw)), kind=4)
                          j16 = 2*(igb + ngb*(ittp-1)) - 1
                          bh(j16)   = real(real(v16), kind=2)
                          bh(j16+1) = real(aimag(v16), kind=2)
                        enddo
                      enddo
                      ierr = cmm_h16_d(wc, bh, fbw, czwc_iw, ngb, ntw, ngb, opa=m_op_C, key = 100000 + iw)
                    else
#endif
                    !$acc parallel loop collapse(2) present(zmel, wgtiw, itw, itpw) async(1)
                    do ittp = 1, ntw
                      do igb = 1, ngb
                        wz_iw(igb,ittp) = cmplx(wgtiw(ittp,iw)*zmel(igb,itw(ittp,iw),itpw(ittp,iw)), kind=kp)
                      enddo
                    enddo
                    ierr = gemm(wc, wz_iw, czwc_iw, ngb, nttp(iw), ngb, opA=m_op_C, key = 100000 + iw, & ! W(omega)
                                policy = BACKEND_SIGMA)
#if defined(__MP) && defined(__GPU)
                    endif
#endif
                    !$acc parallel loop collapse(2) present(itw, itpw) async(1)
                    do ittp = 1, ntw
                      do igb = 1, ngb       ! each (it,itp) appears once per mesh point iw: no two ittp write the same column
                        czmelwc(igb,itw(ittp,iw),itpw(ittp,iw)) = czmelwc(igb,itw(ittp,iw),itpw(ittp,iw)) + czwc_iw(igb,ittp)
                      enddo
                    enddo
                  enddo iwreal
1113              continue !endif
                EndBlock CorrelationSelfEnergyRealAxis
                if (debug) call writemem('    endof CorrelationSelfEnergyRealAxis')
                call stopwatch_pause(sxs_cr)

                !$acc host_data use_device(zmel, zsec)
                ierr = gemm(czmelwc, zmel, zsec, sxs_ntqxx, sxs_ntqxx, nbb*(ns2-ns1+1), opA = m_op_C, beta = CONE, ldC = ntq, &
                            policy = BACKEND_SIGMA)
                !$acc end host_data
                !$acc kernels loop independent async(1)
                do itp = 1, sxs_ntqxx
                  ! zsec(itp,itp) = real(zsec(itp,itp),kind=kp)+img*min(-real((img*zsec(itp,itp)),kind=kp),0_kp) !enforce Imzsec<0 !does not work in intel
                  zsec_img = -real((img*zsec(itp,itp)), kind=kp)
                  if (zsec_img > 0_kp) zsec_img = 0_kp
                  zsec(itp,itp) = real(zsec(itp,itp), kind=kp) + img*zsec_img
                enddo
                !$acc end kernels
                !$acc wait(1)
                call sigma_stream_end()
                !$acc exit data delete(wgtim)
                if (allocated(itw)) then
                  !$acc exit data delete(wgtiw, nttp, itw, itpw)
                  deallocate(itw, itpw, wgtiw)
                endif
                if (allocated(wz_iw)) deallocate(wz_iw, czwc_iw)
                deallocate(wv, wc, czmelwc, wzmel, wgtim)
#if defined(__MP) && defined(__GPU)
                deallocate(bh, wmaxiw)
#endif
                if (ipr) call writemem('    endof CorrelationSelfEnergy')
1114            continue
              endblock get_correlation_block  !end subroutine get_correlation
            endassociate

            call stopwatch_pause(sxs_xc)
            write(stdo,ftox) '    End of icount:', icount ,' of', ncount, &
                 'zmel:', ftof(stopwatch_lap_time(sxs_zmel),4),     '(sec)', &
                 'ec(iaxis,host):', ftof(stopwatch_lap_time(sxs_ci),4),  '(sec)', &
                 'ec(raxis,host):', ftof(stopwatch_lap_time(sxs_cr),4),  '(sec)', &
                 'ec:', ftof(stopwatch_lap_time(sxs_xc),4),         '(sec)', &
                 'setwv:', ftof(stopwatch_elapsed_time(sxs_setwv),4), '(sec)'
                 ! '# of computed real sxs_omega bin:', n_nttp
            call flush(stdo)
          enddo NMBATCHloop
        enddo isploopexternal
      enddo iploopexternal
    enddo irotloop
    ReleaseWV: block !subroutine releasewv()
      if (any(kx == kxc(:))) then
        ! wvi/wvr_upper: device allocatable → deallocate frees GPU memory directly
        if (allocated(wvi_upper)) deallocate(wvi_upper)
        if (allocated(wvr_upper)) deallocate(wvr_upper)
        if (allocated(idx_i)) then
          !$acc exit data delete(idx_i)
          deallocate(idx_i)
        endif
        if (allocated(idx_j)) then
          !$acc exit data delete(idx_j)
          deallocate(idx_j)
        endif
        call wv_close_iq_for_read()
      endif
    end block ReleaseWV !  end subroutine releasewv
  end subroutine sxcf_correlation_step_kx

  real(4) function pow2scale(bound)
    !> 2^(14 - exponent(bound)): a value up to bound, times this, is below 2^14 in FP16 (1 for bound 0).
    real(8), intent(in) :: bound
    pow2scale = 1.0
    if (bound > 0d0) pow2scale = scale(1.0, 14 - exponent(bound))
  end function pow2scale
  subroutine sigma_stream_begin()
    !> One batch of Sigma_c on OpenACC queue 1: its kernels are async(1) and the device products of m_blas (cuBLAS,
    !> realsgemm, GEMMul8) go to the stream of queue 1 too, so they stay in order without a host wait per kernel.
    !> That stream is non-blocking (it does not wait for work on the default stream by itself), so the device is
    !> synchronized first: build_zmel ran on the default stream.  sigma_stream_end after !$acc wait(1).
#ifdef __GPU
    use openacc, only: acc_get_cuda_stream
    use cudafor, only: cudaDeviceSynchronize
    use m_blas, only: cublas_set_stream
    integer :: istat
    istat = cudaDeviceSynchronize()
    call cublas_set_stream(acc_get_cuda_stream(1))
#endif
  end subroutine sigma_stream_begin
  subroutine sigma_stream_end()
#ifdef __GPU
    use cudafor, only: cuda_stream_kind
    use m_blas, only: cublas_set_stream
    call cublas_set_stream(0_cuda_stream_kind)
#endif
  end subroutine sigma_stream_end

  ! Post-kxloop teardown: copy zsecall back to host, deallocate workspace, show timers.
  subroutine sxcf_correlation_finalize()
    !$acc exit data copyout(zsecall)
    !$acc exit data delete(sxs_ekc, sxs_omega, freqx, wx, expa_)
    deallocate(sxs_ekc, sxs_eq, sxs_omega)
    call stopwatch_show(sxs_zmel)
    call stopwatch_show(sxs_ci)
    call stopwatch_show(sxs_cr)
    call stopwatch_show(sxs_pole)
    call stopwatch_show(sxs_xc)
  end subroutine sxcf_correlation_finalize

end module m_sxcf_sc

