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
  use m_blas, only: m_op_c, m_op_n, m_op_t
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
  ! Step WB.3b: working state lives at module scope (ecalj-style singleton).
  ! Lifetime: allocated in init / start of exchange, used by step_kx /
  ! kxloop body, deallocated in finalize / end of exchange.
  ! Stopwatches.
  type(stopwatch) :: sxs_zmel, sxs_xc, sxs_cr, sxs_ci, sxs_setwv
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
  logical,external :: cmdopt0 !we need external here
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
                       alpha = cmplx(-sxs_wkkr,0_kp,kind=kp), beta = CONE, ldC = ntq)
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
    call gwinput_init()
    if (gwinput_loaded) then
       sxs_keepwv = tg_KeepWV
    else
       call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
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
    call stopwatch_init(sxs_cr,    'ec realaxis integral')
    call stopwatch_init(sxs_ci,    'ec imagaxis integral')
    call stopwatch_init(sxs_setwv, 'read wv')
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
    integer, intent(in) :: kx, nspinmx
    real(8), intent(in) :: ef, esmr
    integer :: icount, ns1, ns2, kr, nwxi, ns2r, nwx, izz, n_nttp, tri_idx
    integer :: irot, ip, isp
    real(8) :: q(3), qibz_k(3), qbz_kr(3), qk(3)
    logical :: debug
    real(8), parameter :: ddw = 10d0
    integer, allocatable :: idx_i(:), idx_j(:)
    character(64) :: charli
    character(8)  :: charext
    debug = cmdopt0('--debug')
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
          if (ipr) write(stdo,ftox) 'save WVI and WVR on CPU and GPU (if GPU is used) memory. This requires sufficient memory'
          ! (MO) wvi & wvr are also allocated in CPU memory and are note needed for GPU calcualtion. but allocation of huge
          ! device memory made a error (I don't know the reason). therefore, we used openacc data copyin procedure
          ! but it is usually ok becuase CPU memoery size is always larger than that of GPU.
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
          allocate(wv(nblochpmx,nblochpmx))
          allocate(wvi_upper(ngb*(ngb+1)/2, sxs_wi_ini:sxs_wi_fin))
          allocate(wvr_upper(ngb*(ngb+1)/2, sxs_wr_ini:sxs_wr_fin))
          do iw = sxs_wi_ini, sxs_wi_fin
            if (iw == 0) then
              call wv_get_real(iw, wv)
            else
              call wv_get_imag(iw, wv)
            endif
            do tri_idx = 1, ngb*(ngb+1)/2
              i = idx_i(tri_idx)
              j = idx_j(tri_idx)
              wvi_upper(tri_idx,iw) = wv(i,j)
            enddo
          enddo
          do iw = sxs_wr_ini, sxs_wr_fin
            call wv_get_real(iw, wv)
            do tri_idx = 1, ngb*(ngb+1)/2
              i = idx_i(tri_idx)
              j = idx_j(tri_idx)
              wvr_upper(tri_idx,iw) = (wv(i,j) + conjg(wv(j,i)))*0.5_kp
            enddo
          enddo
          !$acc enter data copyin(wvi_upper, wvr_upper, idx_i, idx_j)
          call stopwatch_pause(sxs_setwv)
          if (ipr) write(stdo, '(X,A,2F8.3)') 'WVI/WVR : sizes (GB)', dble(size(wvi_upper))*kp*2/gb, dble(size(wvr_upper))*kp*2/gb
          call stopwatch_show(sxs_setwv)
          deallocate(wv)
        endif
      endif
      !end subroutine setwv
    end block SetWVblock
    izz = 0
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
          sxs_nt0p = count(sxs_ekc < ef + ddw*esmr)
          sxs_nt0m = count(sxs_ekc < ef - ddw*esmr)
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
                real(8), parameter :: wfaccut=1d-8
                complex(kind=kp), parameter :: img=(0_kp,1_kp)
                complex(kind=kp) :: beta
                complex(kind=kp), allocatable :: czmelwc(:,:,:)
                integer :: it, itp, iw, ierr, i, j
                complex(kind=kp), allocatable :: wv(:,:), wc(:,:)
                real(kind=kp) :: zsec_img
#ifdef __GPU
                attributes(device) :: czmelwc, wc
#endif
                if (ns1 > ns2) goto 1114 !instead of return
                allocate(wv(nblochpmx,nblochpmx))
                allocate(wc(ngb,ngb))
                allocate(czmelwc, mold = zmel)
                call stopwatch_start(sxs_ci)
                CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)! Integration along ImAxis for zwz(sxs_omega)
                  use m_readfreq_r, only: wt=>wwx, x=>freqx
                  real(8):: wgtim_(0:npm*niw), wgtim(0:npm*niw,ns1:ns2,sxs_ntqxx), we, cons(niw), omd(niw), omd2w(niw)
                  real(8):: sig, sig2, aw, aw2
                  integer :: igb
                  complex(kind=kp), allocatable :: wzmel(:,:,:)
#ifdef __GPU
                  attributes(device) :: wzmel
#endif
                  sig = .5d0*esmr
                  sig2 = 2d0*(.5d0*esmr)**2
                  itpdo: do itp = 1, sxs_ntqxx
                   itpo: do it = ns1, ns2
                      we = .5d0*(sxs_omega(itp) - sxs_ekc(it)) !we in hartree unit (atomic unit)
                      aw = abs(ua_*we)
                      aw2 = aw*aw
                      if (it<=nctot) then ! if w = e the integral = -v(0)/2 ! frequency integral
                        cons = 1d0/(we**2*x**2 + (1d0-x)**2)
                        wgtim_(1:niw)= we*cons*wt*(-1d0/pi)
                        wgtim_(0)=merge(-sum(wgtim_(1:niw)*expa_)-0.5d0*dsign(1d0,we)*dexp(we**2*ua_**2)*erfc(ua_*dabs(we)),0d0,&
                             mask=dabs(we)<rmax/ua_)
                        if (npm==2) wgtim_(niw+1:2*niw) = cons*(1d0/x-1d0)*wt/pi !Asymmetric contribution need check
                      else
                        omd   = 1d0/x - 1d0
                        omd2w = omd**2 + we**2
                        where(omd2w/sig2  > 5d-3) cons = (1d0 - exp(-omd2w/sig2))/omd2w
                        where(omd2w/sig2 <= 5d-3) cons = (1d0/sig2 -omd2w/sig2**2/2d0 +omd2w**2/sig2**3/6d0 -omd2w**3/sig2**4/24d0&
                             + omd2w**4/sig2**5/120d0 - omd2w**5/sig2**6/720d0)
                        wgtim_(1:niw) = -we*cons*wt/(x**2)/pi
                        wgtim_(0) =  we*sum(cons*expa_*wt/(x**2))/pi &
                             + dsign(1d0,we)*.5d0*exp(aw2)*( erfc(sqrt(aw2 + we**2/sig2)) - erfc(aw) ) !See Eq.(57) in PRB165106
                        if (npm==2) wgtim_(niw+1:2*niw) = cons*omd*wt/(x**2)/pi !Asymmetric contribution need check
                      endif
                      wgtim(:,it,itp) = sxs_wkkr*wgtim_ !! Integration weight wgtim along im axis for zwz(0:niw*npm)
                    enddo itpo
                  enddo itpdo
                  allocate(wzmel(1:ngb,ns1:ns2,1:sxs_ntqxx))
                  if (debug) call writemem('    Goto iwimag')
                  if (debug) write(stdo,ftox) 'mmmmSc size of mm in imagaxis', (ns2-ns1+1)*sxs_ntqxx, ngb, ngb
                  !$acc data copyin(wgtim)
                  iwimag: do iw = sxs_wi_ini, sxs_wi_fin ! iwimag:do iw = 0, niw !niw is ~10. ixx=0 is for sxs_omega=0 nw_i=0 (Time reversal) or nw_i =-nw
                    if (iw < 0 .or. iw > niw) cycle
                    call stopwatch_start(sxs_setwv)
                    if (sxs_keepwv) then
                      !$acc kernels loop independent present(wvi_upper, idx_i, idx_j)
                      do tri_idx = 1, ngb*(ngb+1)/2
                        i = idx_i(tri_idx)
                        j = idx_j(tri_idx)
                        wc(i,j) = wvi_upper(tri_idx,iw)
                        wc(j,i) = conjg(wvi_upper(tri_idx,iw))
                      enddo
                      !$acc end kernels
                    else
                      if (iw == 0) call wv_get_real(iw, wv)
                      if (iw > 0)  call wv_get_imag(iw, wv)
                      wc(1:ngb,1:ngb) = wv(1:ngb,1:ngb)  !copy to GPU
                    endif
                    call stopwatch_pause(sxs_setwv)
                    !$acc kernels loop independent collapse(2) present(zmel)
                    do itp = 1, sxs_ntqxx
                      do it = ns1, ns2
                        wzmel(1:ngb,it,itp) = cmplx(wgtim(iw,it,itp)*zmel(1:ngb,it,itp), kind=kp)
                      enddo
                    enddo
                    !$acc end kernels
                    !the most time-consuming part in the correlation part
                    beta = CONE
                    if (iw == sxs_wi_ini) beta = CZERO
                    ierr = gemm(wc, wzmel, czmelwc, ngb, (ns2-ns1+1)*sxs_ntqxx, ngb, beta = beta, opA = m_op_C)
                  enddo iwimag
                  !$acc end data
                  deallocate(wzmel)
                EndBlock CorrelationSelfEnergyImagAxis
                if (debug) call writemem('    endof CorrelationSelfEnergyImagAxis')
                call stopwatch_pause(sxs_ci)

                call stopwatch_start(sxs_cr)
                CorrelationSelfEnergyRealAxis: Block !Real Axis integral. Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
                  use m_wfac, only: wfacx2, weavx2
                  integer :: itini, itend, ittp, ittp3(3), ixs, nttp_max, nttp(0:nw), i, j
                  real(8) :: we_(ns1:ns2r,sxs_ntqxx), wfac_(ns1:ns2r,sxs_ntqxx), omg, amat(3,3), wgt3ititp(3)
                  complex(kind=kp), allocatable :: wz_iw(:,:), czwc_iw(:,:)
                  real(8), allocatable :: wgtiw(:,:)
                  integer, allocatable :: itw(:,:), itpw(:,:)
#ifdef __GPU
                  attributes(device) :: wz_iw, czwc_iw
#endif
                  nttp = 0
                  itploop: do itp = 1, sxs_ntqxx
                    omg   = sxs_omega(itp)
                    itini = merge(max(ns1,sxs_nt0m+1),  ns1, mask= omg>=ef)
                    itend = merge(ns2r,  min(sxs_nt0p,ns2r), mask= omg>=ef)
                    do it = itini, itend
                      wfac_(it,itp) = wfacx2(omg, ef, sxs_ekc(it), merge(0d0,esmr,mask=it<=nctot))
                      if (wfac_(it,itp) < wfaccut) cycle
                      we_(it,itp)  = .5d0*abs(omg - weavx2(omg,ef, sxs_ekc(it),esmr))
                      ixs = findloc(freq_r(1:nw)>we_(it,itp), value=.true., dim=1)
                      nttp(ixs-1:ixs+1) = nttp(ixs-1:ixs+1) + 1
                    enddo
                  enddo itploop
                  nttp_max = maxval(nttp)
                  if (nttp_max <= 0) goto 1113
                  allocate (itw(nttp_max,0:nw),  source = 0)
                  allocate (itpw(nttp_max,0:nw), source = 0)
                  allocate (wgtiw(nttp_max,0:nw),source = 0d0)
                  nttp = 0
                  itploopFORwgtiw: do itp = 1, sxs_ntqxx
                    omg   = sxs_omega(itp)
                    itini = merge(max(ns1,sxs_nt0m+1),  ns1, mask= omg>=ef)
                    itend = merge(ns2r,  min(sxs_nt0p,ns2r), mask= omg>=ef)
                    do it = itini, itend     ! sxs_nt0p corresponds to efp
                      wfac_(it,itp) = wfacx2(omg, ef, sxs_ekc(it), merge(0d0,esmr,mask=it<=nctot)) !Gaussian smearing
                      if (wfac_(it,itp) < wfaccut) cycle
                      wfac_(it,itp) =  wfac_(it,itp)*sxs_wkkr*dsign(1d0, omg-ef) !wfac_ = $w$ weight (smeared thus truncated by ef). See the sentences.
                      we_(it,itp)   = .5d0*abs(omg - weavx2(omg,ef, sxs_ekc(it),esmr)) !we_= \bar{\omega_\epsilon} in sentences next to Eq.58 in PRB76,165106 (2007)
                      ixs = findloc(freq_r(1:nw)>we_(it,itp), value=.true., dim=1)
                      associate(x => we_(it,itp), xi => freq_r(ixs-1:ixs+1)) !x=>we_ is \omega_\epsilon in Eq.(55).
                        amat(1:3,1) = 1d0                 !old version: call alagr3z2wgt(we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
                        amat(1:3,2) = xi(1:3)**2
                        amat(1:3,3) = xi(1:3)**4
                        wgt3ititp = wfac_(it,itp)*matmul([1d0, x**2, x**4], inverse33(amat))
                      end associate
                      nttp(ixs-1:ixs+1) = nttp(ixs-1:ixs+1) + 1
                      ittp3(1:3) = nttp(ixs-1:ixs+1)
                      forall(i=1:3) itw(ittp3(i),  ixs-2+i) = it
                      forall(i=1:3) itpw(ittp3(i), ixs-2+i) = itp
                      forall(i=1:3) wgtiw(ittp3(i),ixs-2+i) = wgt3ititp(i)
                    enddo
                  enddo itploopFORwgtiw
                  n_nttp = count(nttp(sxs_wr_ini:sxs_wr_fin) > 0)
                  allocate(wz_iw(ngb,nttp_max), czwc_iw(ngb,nttp_max))
                  !$acc data copyin(wgtiw, nttp, itw, itpw)
                  iwreal: do iw = sxs_wr_ini, sxs_wr_fin
                    if (iw < nwxi .or. iw > nwx) cycle
                    if (nttp(iw) < 1) cycle
                    call stopwatch_start(sxs_setwv)
                    if (sxs_keepwv) then
                      !$acc kernels loop independent present(wvr_upper, idx_i, idx_j)
                      do tri_idx = 1, ngb*(ngb+1)/2
                        i = idx_i(tri_idx)
                        j = idx_j(tri_idx)
                        wc(i,j) = wvr_upper(tri_idx,iw)
                        wc(j,i) = conjg(wvr_upper(tri_idx,iw))
                      enddo
                      !$acc end kernels
                    else
                      call wv_get_real(iw, wv)
                      wc(1:ngb,1:ngb) = wv(1:ngb,1:ngb)  !copy to GPU
                      !$acc kernels
                      wc(:,:) = (wc(:,:) + transpose(conjg(wc(:,:))))*0.5_kp
                      !$acc end kernels
                    endif
                    call stopwatch_pause(sxs_setwv)
                    !$acc kernels loop independent present(zmel)
                    do ittp = 1, nttp(iw)
                      it = itw(ittp,iw); itp = itpw(ittp,iw)
                      wz_iw(1:ngb,ittp) = cmplx(wgtiw(ittp,iw)*zmel(1:ngb,it,itp), kind=kp)
                    enddo
                    !$acc end kernels
                    ierr = gemm(wc, wz_iw, czwc_iw, ngb, nttp(iw), ngb, opA=m_op_C)
                    !$acc kernels loop independent
                    do ittp = 1, nttp(iw)
                      it = itw(ittp,iw); itp = itpw(ittp,iw)
                      czmelwc(1:ngb,it,itp) = czmelwc(1:ngb,it,itp) + czwc_iw(1:ngb,ittp)
                    enddo
                    !$acc end kernels
                  enddo iwreal
                  !$acc end data
                  deallocate(wz_iw, czwc_iw)
1113              continue !endif
                EndBlock CorrelationSelfEnergyRealAxis
                if (debug) call writemem('    endof CorrelationSelfEnergyRealAxis')
                call stopwatch_pause(sxs_cr)

                !$acc host_data use_device(zmel, zsec)
                ierr = gemm(czmelwc, zmel, zsec, sxs_ntqxx, sxs_ntqxx, nbb*(ns2-ns1+1), opA = m_op_C, beta = CONE, ldC = ntq)
                !$acc end host_data
                !$acc kernels loop independent
                do itp = 1, sxs_ntqxx
                  ! zsec(itp,itp) = real(zsec(itp,itp),kind=kp)+img*min(-real((img*zsec(itp,itp)),kind=kp),0_kp) !enforce Imzsec<0 !does not work in intel
                  zsec_img = -real((img*zsec(itp,itp)), kind=kp)
                  if (zsec_img > 0_kp) zsec_img = 0_kp
                  zsec(itp,itp) = real(zsec(itp,itp), kind=kp) + img*zsec_img
                enddo
                !$acc end kernels
                deallocate(wv, wc, czmelwc)
                if (ipr) call writemem('    endof CorrelationSelfEnergy')
1114            continue
              endblock get_correlation_block  !end subroutine get_correlation
            endassociate

            call stopwatch_pause(sxs_xc)
            write(stdo,ftox) '    End of icount:', icount ,' of', ncount, &
                 'zmel:', ftof(stopwatch_lap_time(sxs_zmel),4),     '(sec)', &
                 'ec(iaxis):', ftof(stopwatch_lap_time(sxs_ci),4),  '(sec)', &
                 'ec(raxis):', ftof(stopwatch_lap_time(sxs_cr),4),  '(sec)', &
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
        if (allocated(wvi_upper)) then
          !$acc exit data delete(wvi_upper)
          deallocate(wvi_upper)
        endif
        if (allocated(wvr_upper)) then
          !$acc exit data delete(wvr_upper)
          deallocate(wvr_upper)
        endif
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

  ! Post-kxloop teardown: copy zsecall back to host, deallocate workspace, show timers.
  subroutine sxcf_correlation_finalize()
    !$acc exit data copyout(zsecall)
    deallocate(sxs_ekc, sxs_eq, sxs_omega)
    call stopwatch_show(sxs_zmel)
    call stopwatch_show(sxs_ci)
    call stopwatch_show(sxs_cr)
    call stopwatch_show(sxs_xc)
  end subroutine sxcf_correlation_finalize

  pure function inverse33(matrix) result(inverse) !Inverse of 3X3 matrix
    implicit none
    real(8),intent(in) :: matrix(3,3)
    real(8) :: inverse(3,3), det
    inverse(:,1)= crossf(matrix(:,2),matrix(:,3))
    inverse(:,2)= crossf(matrix(:,3),matrix(:,1))
    inverse(:,3)= crossf(matrix(:,1),matrix(:,2))
    det = sum(matrix(:,1)*inverse(:,1))
    inverse = transpose(inverse)
    inverse = inverse/det
  end function inverse33
  pure function crossf(a,b) result(c)
    implicit none
    intent(in):: a,b
    real(8):: a(3),b(3),c(3)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
  end function crossf
end module m_sxcf_sc
