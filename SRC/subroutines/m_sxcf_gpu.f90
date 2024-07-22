module m_sxcf_gpu
  use m_readeigen, only: Readeval
  use m_zmel, only: get_zmel_init => get_zmel_init_gpu, Setppovlz, zmel, nbb
  use m_itq, only: itq, ntq, nbandmx
  use m_genallcf_v3, only: nlmto, nspin, nctot, niw, ecore !,symgg
  use m_read_bzdata, only: qibz, qbz, wk=>wbz, nqibz, nqbz, wklm, lxklm, wqt=>wt
  use m_readVcoud, only: Readvcoud, vcoud, ngb, ngc
  use m_readfreq_r, only: freq_r, nw_i, nw, freqx, wx=>wwx, nblochpmx, mrecl, expa_, npm
  use m_readhbe, only: nband
  use m_readgwinput, only: ua_, corehole, wcorehole
  use m_ftox
  use m_lgunit, only: stdo
  use m_sxcf_count,only: kxc, nstti, nstte, nstte2, nwxic, nwxc, icountini, icountend, irkip, ncount
  use m_nvfortran, only: findloc
  use m_hamindex, only: ngrp
  use m_blas, only: m_op_c, m_op_n, m_op_t
  use m_blas, only: gemm => zmm, gemm_batch => zmm_batch
  use m_kind, only: kp => kindgw
  use m_sxcf_main,only: zsecall
  use m_stopwatch
  use m_mem,only:writemem
  implicit none
  public sxcf_scz_correlation, sxcf_scz_exchange
  private
  real(8), parameter :: pi = 4d0*datan(1d0), fpi = 4d0*pi
  logical, parameter :: timemix = .true.
  complex(kind=kp), parameter :: CONE = (1_kp, 0_kp), CZERO = (0_kp, 0_kp)
  integer :: kx, irot, ip, isp, ntqxx, nt0p, nt0m, ifrcw, ifrcwi 
  real(8) :: wkkr
  integer, allocatable :: ndiv(:), nstatei(:,:), nstatee(:,:)
  real(8), allocatable :: ekc(:), eq(:), omega(:)
  logical :: emptyrun, keepwv
  complex(kind=kp), allocatable :: wvr(:,:,:), wvi(:,:,:)
#ifdef __GPU
  attributes(device) :: wvr, wvi
#endif
  logical, external :: cmdopt0
  type(stopwatch) :: t_sw_zmel, t_sw_xc, t_sw_cr, t_sw_ci
contains
  subroutine sxcf_scz_exchange(ef, esmr, ixc, nspinmx) !ixc is dummy
    implicit none
    integer :: icount, ns1, ns2, kr,izz
    logical, parameter :: debug=.false.
    integer, intent(in) :: nspinmx, ixc
    real(8), intent(in) :: ef, esmr
    real(8) :: q(3), qibz_k(3), qbz_kr(3), qk(3)
    character(64):: charli
    character(8):: charext
    allocate(ekc(nctot+nband), eq(nband)) 
    emptyrun = cmdopt0('--emptyrun')
    if(nw_i/=0) call rx('Current version we assume nw_i=0. Time-reversal symmetry')
    allocate(zsecall(ntq,ntq,nqibz,nspinmx),source=(0d0,0d0)) 
    LoopScheduleCheck: block
    izz=0
    kxloopX:                   do kx  =1,nqibz  
       irotloopX:              do irot=1,ngrp    
          iploopexternalX:     do ip=1,nqibz     
             isploopexternalX: do isp=1,nspinmx  
                kr = irkip(isp,kx,irot,ip)
                if(kr==0) cycle
                NMBATCHloopX: do icount = icountini(isp,ip,irot,kx),icountend(isp,ip,irot,kx) !batch of middle states.
                   izz=izz+1
                   write(stdo,ftox)'- kxloop Schedule ',izz, ' iqibz irot ip isp icount=',kx,irot,ip,isp,icount
                enddo NMBATCHloopX
             enddo isploopexternalX
          enddo iploopexternalX
       enddo irotloopX
    enddo kxloopX
    end block LoopScheduleCheck

    call stopwatch_init(t_sw_zmel, 'zmel')
    call stopwatch_init(t_sw_xc, 'ex')
    !$acc data copyout(zsecall)
    !$acc kernels
      zsecall(1:ntq,1:ntq,1:nqibz,1:nspinmx) = CZERO
    !$acc end kernels
    kxloop: do kx=1, nqibz                         ! kx is irreducible !kx is main axis where we calculate W(kx).
      qibz_k = qibz(:,kx)
      call Readvcoud(qibz_k, kx, NoVcou=.false.)   !Readin ngc,ngb,vcoud ! Coulomb matrix
      call Setppovlz(qibz_k, matz=.true., npr=ngb) !Set ppovlz overlap matrix used in Get_zmel_init in m_zmel
      irotloop: do irot=1, ngrp                    ! (kx,irot) determines qbz(:,kr), which is in FBZ. W(kx) is rotated to be W(g(kx))
        iploopexternal: do ip=1, nqibz             !external index for q of \Sigma(q,isp)
          isploopexternal: do isp=1, nspinmx       !external index
            kr = irkip(isp,kx,irot,ip)
            if(kr==0) cycle
            q = qibz(:,ip)
            qbz_kr = qbz (:,kr)   !rotated qbz vector.
            qk = q - qbz_kr       !<M(qbz_kr) phi(q-qbz_kr)|phi(q)>
            eq = readeval(q,isp)  !readin eigenvalue
            ekc(1:nctot+nband) = [ecore(1:nctot,isp),readeval(qk, isp)]
            ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
            wkkr = wk(kr)
            NMBATCHloop: do icount = icountini(isp,ip,irot,kx), icountend(isp,ip,irot,kx) !batch of middle states.
              ns1 = nstti(icount)  !Range of middle states is [ns1:ns2] for given icount
              ns2 = nstte(icount)  ! 
              call stopwatch_start(t_sw_zmel)
              izz=izz+1
              call writemem(' -- KXloop '//trim(charext(izz))//' iqiqz irot ip isp icount= '//&
                   trim(charli([kx,irot,ip,isp,icount],5)))
              call get_zmel_init(q, qibz_k, irot, qbz_kr, ns1, ns2, isp, 1, ntqxx, isp, nctot, ncc=0, iprx=debug, zmelconjg=.false.)
              call writemem(' --        end of zmel')
              call stopwatch_pause(t_sw_zmel)
              call stopwatch_start(t_sw_xc)
              call get_exchange(ef, esmr, ns1, ns2, zsecall(1,1,ip,isp))
              call stopwatch_pause(t_sw_xc)
              write(stdo,ftox) 'end of icount:', icount ,' of', ncount, &
                               'zmel:', ftof(stopwatch_lap_time(t_sw_zmel),4), '(sec)', &
                               'exch:', ftof(stopwatch_lap_time(t_sw_xc),4),   '(sec)'
              call flush(stdo)
              call writemem('endof ExchangeSelfEnergy')
            enddo NMBATCHloop
          enddo isploopexternal
        enddo iploopexternal
      enddo irotloop
    enddo kxloop
    !$acc end data
    deallocate(ekc, eq)
    call stopwatch_show(t_sw_zmel)
    call stopwatch_show(t_sw_xc)
  endsubroutine sxcf_scz_exchange

  subroutine sxcf_scz_correlation(ef, esmr, ixc, nspinmx) 
    implicit none
    integer, intent(in) :: nspinmx, ixc
    real(8), intent(in) :: ef, esmr
    integer :: icount, ns1, ns2, kr, nwxi, nws, ns2r, nwx,izz
    real(8) :: q(3), qibz_k(3), qbz_kr(3), qk(3)
    logical, parameter :: debug=.false.
    real(8),parameter :: ddw=10d0
    character(64):: charli
    character(8):: charext
    allocate(ekc(nctot+nband), eq(nband), omega(ntq)) 
    emptyrun = cmdopt0('--emptyrun')
    keepwv = cmdopt0('--keepwv')
    if(nw_i/=0) call rx('Current version we assume nw_i=0. Time-reversal symmetry')
    LoopScheduleCheck: block
    izz=0
    kxloopX:                   do kx  =1,nqibz  
       irotloopX:              do irot=1,ngrp    
          iploopexternalX:     do ip=1,nqibz     
             isploopexternalX: do isp=1,nspinmx  
                kr = irkip(isp,kx,irot,ip)
                if(kr==0) cycle
                NMBATCHloopX: do icount = icountini(isp,ip,irot,kx),icountend(isp,ip,irot,kx) !batch of middle states.
                   izz=izz+1
                   write(stdo,ftox)'- KXloop Scheduling ',izz,' iqiqz irot ip isp icount=',kx,irot,ip,isp,icount
                enddo NMBATCHloopX
             enddo isploopexternalX
          enddo iploopexternalX
       enddo irotloopX
    enddo kxloopX
    end block LoopScheduleCheck
    call stopwatch_init(t_sw_zmel, 'zmel')
    call stopwatch_init(t_sw_xc, 'ec')
    call stopwatch_init(t_sw_cr, 'ec realaxis integral')
    call stopwatch_init(t_sw_ci, 'ec imagaxis integral')
    allocate(zsecall(ntq,ntq,nqibz,nspinmx),source=(0d0,0d0)) 
    !$acc data copyout(zsecall)
    !$acc kernels
      zsecall(1:ntq,1:ntq,1:nqibz,1:nspinmx) = CZERO
    !$acc end kernels
    kxloop: do kx=1, nqibz                         ! kx is irreducible !kx is main axis where we calculate W(kx).
      qibz_k = qibz(:,kx)
      call Readvcoud(qibz_k, kx, NoVcou=.false.)   !Readin ngc,ngb,vcoud ! Coulomb matrix
      call Setppovlz(qibz_k, matz=.true., npr=ngb) !Set ppovlz overlap matrix used in Get_zmel_init in m_zmel
      call setwv()
      irotloop: do irot=1, ngrp                    ! (kx,irot) determines qbz(:,kr), which is in FBZ. W(kx) is rotated to be W(g(kx))
        iploopexternal: do ip=1, nqibz             !external index for q of \Sigma(q,isp)
          isploopexternal: do isp=1, nspinmx       !external index
            kr = irkip(isp,kx,irot,ip)
            if(kr==0) cycle
            q = qibz(:,ip)
            qbz_kr = qbz (:,kr)   !rotated qbz vector.
            qk = q - qbz_kr       !<M(qbz_kr) phi(q-qbz_kr)|phi(q)>
            eq = readeval(q,isp)  !readin eigenvalue
            wkkr = wk(kr)
            ekc(1:nctot+nband) = [ecore(1:nctot,isp),readeval(qk, isp)]
            ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
            omega(1:ntq) = eq(1:ntq)  
            nt0p = count(ekc<ef+ddw*esmr)  
            nt0m = count(ekc<ef-ddw*esmr) 
            NMBATCHloop: do icount = icountini(isp,ip,irot,kx), icountend(isp,ip,irot,kx) !batch of middle states.
              ns1 = nstti(icount)  !Range of middle states is [ns1:ns2] for given icount
              ns2 = nstte(icount)  ! 
              nwxi = nwxic(icount)  !minimum omega for W
              nwx = nwxc(icount)   !max omega for W
              ns2r = nstte2(icount) !Range of middle states [ns1:ns2r] for CorrelationSelfEnergyRealAxis
              call writemem(' -- KXloop '//trim(charext(izz))//' iqiqz irot ip isp icount= '//&
                   trim(charli([kx,irot,ip,isp,icount],5)))
              call stopwatch_start(t_sw_zmel)
              call get_zmel_init(q, qibz_k, irot, qbz_kr, ns1, ns2, isp, 1, ntqxx, isp, nctot, ncc=0, iprx=debug, zmelconjg=.false.)
              call stopwatch_pause(t_sw_zmel)
              call stopwatch_start(t_sw_xc)
              call get_correlation(ef, esmr, ns1, ns2, ns2r, nwxi, nwx, zsecall(1,1,ip,isp))
              call stopwatch_pause(t_sw_xc)
              write(stdo,ftox) 'end of icount:', icount ,' of', ncount, &
                               'zmel:', ftof(stopwatch_lap_time(t_sw_zmel),4), '(sec)', &
                               'ec(iaxis):', ftof(stopwatch_lap_time(t_sw_ci),4),   '(sec)', &
                               'ec(raxis):', ftof(stopwatch_lap_time(t_sw_cr),4),   '(sec)', &
                               'ec:', ftof(stopwatch_lap_time(t_sw_xc),4),   '(sec)'
              call flush(stdo)
            enddo NMBATCHloop
          enddo isploopexternal
        enddo iploopexternal
      enddo irotloop
      call releasewv()
    enddo kxloop
    !$acc end data
    deallocate(ekc, eq, omega)
    call stopwatch_show(t_sw_zmel)
    call stopwatch_show(t_sw_ci)
    call stopwatch_show(t_sw_cr)
    call stopwatch_show(t_sw_xc)
  endsubroutine sxcf_scz_correlation

  subroutine get_exchange(ef, esmr, ns1, ns2, zsec)
    implicit none
    integer, intent(in) :: ns1, ns2
    real(8), intent(in) :: ef, esmr
    complex(kind=kp), intent(inout) ::zsec(ntq,ntq)
    real(8) :: wfacx ! external function
    real(8), allocatable :: vcoud_buf(:), wtff(:)
    complex(kind=kp), allocatable :: vzmel(:,:,:)
    integer :: it, itp, itpp, ierr
#ifdef __GPU
    attributes(device) :: vzmel, vcoud_buf
#endif
    if(ns1 > ns2) return
    allocate(wtff(ns1:ns2))
    wtff(ns1:nctot) = 1d0
    do it = max(nctot+1,ns1), ns2
      wtff(it) = wfacx(-1d99, ef, ekc(it), esmr)
    enddo
    if(corehole) wtff(ns1:nctot) = wtff(ns1:nctot) * wcorehole(ns1:nctot,isp)
    allocate(vcoud_buf(ngb))
    !$acc data copyin(vcoud, wklm(1), wk(1), wtff, zmel) copy(zsec)
    !$acc kernels
    vcoud_buf(1:ngb) = vcoud(1:ngb)
    !$acc end kernels
    if(kx == 1) vcoud_buf(1) = wklm(1)*fpi*sqrt(fpi)/wk(1)
    allocate(vzmel(1:nbb,ns1:ns2,1:ntqxx))
    !$acc kernels loop independent collapse(2)
    do itpp = 1, ntqxx
      do it = ns1, ns2 
        vzmel(1:ngb,it,itpp) = wtff(it)*vcoud_buf(1:ngb)*zmel(1:ngb,it,itpp)
      enddo
    enddo
    !$acc end kernels
    !$acc host_data use_device(zmel)
    ierr = gemm(zmel, vzmel, zsec, ntqxx, ntqxx, (ns2-ns1+1)*nbb, opA = m_op_c, &
                alpha = cmplx(-wkkr,0_kp,kind=kp), beta = CONE, ldC = ntq)
    !$acc end host_data
    !$acc end data
    deallocate(vzmel, vcoud_buf, wtff)
  end subroutine get_exchange

  subroutine get_correlation(ef, esmr, ns1, ns2, ns2r, nwxi, nwx, zsec)
    integer, intent(in) :: ns1, ns2, nwxi, nwx, ns2r
    real(8), intent(in) :: ef, esmr
    complex(kind=kp), intent(inout) :: zsec(ntq,ntq)
    real(8), parameter :: wfaccut=1d-8
    complex(kind=kp), parameter :: img=(0_kp,1_kp)
    complex(kind=kp) :: beta
    complex(kind=kp), allocatable :: czmelwc(:,:,:)
    integer :: it, itp, iw, ierr, i
    complex(kind=kp), allocatable :: wv(:,:), wc(:,:)
#ifdef __GPU
    attributes(device) :: czmelwc, wc
#endif

    if(ns1 > ns2) return
    allocate(wv(nblochpmx,nblochpmx))
    allocate(wc, mold = wv)
    call stopwatch_start(t_sw_ci)
    CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)! Integration along ImAxis for zwz(omega) 
      use m_readfreq_r, only: wt=>wwx, x=>freqx
      real(8), parameter:: rmax=2d0
      real(8):: wgtim_(0:npm*niw), wgtim(0:npm*niw,ntqxx,ns1:ns2), we, cons(niw), omd(niw), omd2w(niw)
      real(8):: sig, sig2, aw, aw2
      integer :: igb
      complex(kind=kp), allocatable :: wzmel(:,:,:),  czwc(:,:,:)
#ifdef __GPU
      attributes(device) :: wzmel, czwc
#endif
      sig = .5d0*esmr
      sig2 = 2d0*(.5d0*esmr)**2
      do it = ns1, ns2
        do itp = 1, ntqxx
          we = .5d0*(omega(itp)-ekc(it)) !we in hartree unit (atomic unit)
          aw = abs(ua_*we)
          aw2 = aw*aw
          if(it<=nctot) then ! if w = e the integral = -v(0)/2 ! frequency integral
            cons = 1d0/(we**2*x**2 + (1d0-x)**2)
            wgtim_(1:niw)= we*cons*wt*(-1d0/pi)
            wgtim_(0)=merge(-sum(wgtim_(1:niw)*expa_)-0.5d0*dsign(1d0,we)*dexp(we**2*ua_**2)*erfc(ua_*dabs(we)), 0d0,&
                 mask=dabs(we)<rmax/ua_)
            if(npm==2) wgtim_(niw+1:2*niw) = cons*(1d0/x-1d0)*wt/pi !Asymmetric contribution need check
          else
            omd   = 1d0/x - 1d0
            omd2w = omd**2 + we**2
            where(omd2w/sig2  > 5d-3) cons=(1d0 - exp (- omd2w/sig2))/omd2w
            where(omd2w/sig2 <= 5d-3) cons=(1d0/sig2 - omd2w/sig2**2/2d0 +omd2w**2/sig2**3/6d0 -omd2w**3/sig2**4/24d0 &
                 + omd2w**4/sig2**5/120d0- omd2w**5/sig2**6/720d0 )
            wgtim_(1:niw) = -we*cons*wt/(x**2)/pi
            wgtim_(0) =  we*sum(cons*expa_*wt/(x**2))/pi &
                 + dsign(1d0,we)*.5d0*exp(aw2)*( erfc(sqrt(aw2 + we**2/sig2)) -erfc(aw) ) !Gaussian part erfc(2023feb)
            if(npm==2) wgtim_(niw+1:2*niw) = cons*omd*wt/(x**2)/pi !Asymmetric contribution need chack
          endif
          wgtim(:,itp,it)= wkkr*wgtim_ !! Integration weight wgtim along im axis for zwz(0:niw*npm)
        enddo 
      enddo 
      allocate(wzmel(1:ngb,ns1:ns2,1:ntqxx), czwc(ns1:ns2,1:ntqxx,1:ngb))

      call writemem('Goto iwimag')
      !$acc data copyin(wgtim, zmel)
      iwimag:do iw = 0, niw !niw is ~10. ixx=0 is for omega=0 nw_i=0 (Time reversal) or nw_i =-nw
        if(emptyrun) cycle
        if(keepwv) then
          if(iw == 0) wc = wvr(:,:,iw)
          if(iw > 0) wc = wvi(:,:,iw)
        else
          if(iw == 0) read(ifrcw,rec=1+(0-nw_i)) wv!direct access Wc(0) = W(0)-v ! nw_i=0 (Time reversal) or nw_i =-nw
          if(iw > 0) read(ifrcwi,rec=iw) wv ! direct access read Wc(i*omega)=W(i*omega)-v
          wc = wv
        endif
        !$acc kernels loop independent collapse(2)
        do itp = 1, ntqxx
          do it = ns1, ns2
            wzmel(1:ngb,it,itp) = wgtim(iw,itp,it)*zmel(1:ngb,it,itp)
          enddo
        enddo
        !$acc end kernels
        !the most time-consuming part in the correlation part
        beta = CONE
        if(iw == 0) beta = CZERO
        ierr = gemm(wzmel, wc, czwc, (ns2-ns1+1)*ntqxx, ngb, ngb, opA = m_op_c, beta = beta, ldB = nblochpmx)
      enddo iwimag
      !$acc end data
      deallocate(wzmel)
      allocate(czmelwc(1:nbb,ns1:ns2,1:ntqxx)) !same size with zmel
      !$acc kernels loop independent
      do igb = 1, ngb
        czmelwc(igb,ns1:ns2,1:ntqxx) = czwc(ns1:ns2,1:ntqxx,igb)
      enddo
      !$acc end kernels
      deallocate(czwc)
    EndBlock CorrelationSelfEnergyImagAxis
    call writemem('endof CorrelationSelfEnergyImagAxis')
    call stopwatch_pause(t_sw_ci)

    call stopwatch_start(t_sw_cr)
    CorrelationSelfEnergyRealAxis: Block !Real Axis integral. Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
      use m_wfac, only: wfacx2, weavx2
      integer :: itini, itend, ittp, ittp3(3), ixs, nttp_max, nttp(0:nw)
      real(8) :: we_(ns1:ns2r,ntqxx), wfac_(ns1:ns2r,ntqxx), omg, amat(3,3), wgt3ititp(3)
      complex(kind=kp), allocatable :: wz_iw(:,:), czwc_iw(:,:)
      real(8), allocatable :: wgtiw(:,:)
      integer, allocatable :: itw(:,:), itpw(:,:)
#ifdef __GPU 
      attributes(device) :: wz_iw, czwc_iw
#endif
      nttp = 0
      do itp = 1, ntqxx
        omg  = omega(itp)
        itini = merge(max(ns1,nt0m+1),  ns1, mask= omg>=ef)
        itend = merge(ns2r,  min(nt0p,ns2r), mask= omg>=ef)
        do it = itini, itend
          wfac_(it,itp) = wfacx2(omg, ef, ekc(it), merge(0d0,esmr,mask=it<=nctot))
          if(wfac_(it,itp) < wfaccut) cycle 
          we_(it,itp)  = .5d0*abs(omg-weavx2(omg,ef, ekc(it),esmr))
          ixs = findloc(freq_r(1:nw)>we_(it,itp),value=.true.,dim=1)
          nttp(ixs-1:ixs+1) = nttp(ixs-1:ixs+1) + 1
        enddo
      enddo
      nttp_max = maxval(nttp)
      if(nttp_max > 0) then 
        allocate (itw(nttp_max,0:nw), source = 0)
        allocate (itpw(nttp_max,0:nw), source = 0)
        allocate (wgtiw(nttp_max,0:nw), source = real(0,kind=kp))
        nttp = 0
        do itp = 1, ntqxx
          omg  = omega(itp)
          itini = merge(max(ns1,nt0m+1),  ns1, mask= omg>=ef)
          itend = merge(ns2r,  min(nt0p,ns2r), mask= omg>=ef)
          do it = itini, itend     ! nt0p corresponds to efp
            wfac_(it,itp) = wfacx2(omg, ef, ekc(it), merge(0d0,esmr,mask=it<=nctot)) !Gaussian smearing 
            if(wfac_(it,itp)<wfaccut) cycle 
            wfac_(it,itp)=  wfac_(it,itp)*wkkr*dsign(1d0,omg-ef) !wfac_ = $w$ weight (smeared thus truncated by ef). See the sentences.
            we_(it,itp)  = .5d0*abs(omg-weavx2(omg,ef, ekc(it),esmr)) !we_= \bar{\omega_\epsilon} in sentences next to Eq.58 in PRB76,165106 (2007)
            ixs = findloc(freq_r(1:nw)>we_(it,itp),value=.true.,dim=1)
            associate(x=>we_(it,itp),xi=>freq_r(ixs-1:ixs+1))!x=>we_ is \omega_\epsilon in Eq.(55). 
              amat(1:3,1)= 1d0                 !old version: call alagr3z2wgt(we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
              amat(1:3,2)= xi(1:3)**2
              amat(1:3,3)= xi(1:3)**4
              wgt3ititp = wfac_(it,itp)*matmul([1d0,x**2,x**4], inverse33(amat)) 
            endassociate
            nttp(ixs-1:ixs+1) = nttp(ixs-1:ixs+1) + 1
            ittp3(1:3) = nttp(ixs-1:ixs+1)
            forall(i=1:3) itw(ittp3(i), ixs-2+i) = it
            forall(i=1:3) itpw(ittp3(i), ixs-2+i) = itp
            forall(i=1:3) wgtiw(ittp3(i), ixs-2+i) = wgt3ititp(i)
          enddo
        enddo
        allocate(wz_iw(ngb,nttp_max), czwc_iw(nttp_max,ngb))
        !$acc data copyin(wgtiw, nttp, itw, itpw)
        do iw = nwxi, nwx
          if(nttp(iw) < 1) cycle
          if(keepwv) then
            !$acc kernels
            wc(:,:) = (wvr(:,:,iw) + transpose(dconjg(wvr(:,:,iw))))*0.5d0
            !$acc end kernels
          else
            read(ifrcw,rec=iw-nw_i+1) wv
            wc = (wv + transpose(dconjg(wv)))*0.5d0  !copy to GPU
          endif
          !$acc kernels loop independent
          do ittp = 1, nttp(iw) 
            it = itw(ittp,iw); itp = itpw(ittp,iw)
            wz_iw(1:ngb,ittp) = wgtiw(ittp,iw)*zmel(1:ngb,it,itp)
          enddo
          !$acc end kernels
          ierr = gemm(wz_iw, wc, czwc_iw, nttp(iw), ngb, ngb, opA = m_op_c, ldB = nblochpmx, ldC = nttp_max)
          !$acc kernels loop independent
          do ittp = 1, nttp(iw) 
            it = itw(ittp,iw); itp = itpw(ittp,iw)
            czmelwc(1:ngb,it,itp) = czmelwc(1:ngb,it,itp) + czwc_iw(ittp,1:ngb)
          enddo
          !$acc end kernels
        enddo
        !$acc end data
        deallocate(wz_iw, czwc_iw)
      endif
    EndBlock CorrelationSelfEnergyRealAxis
    call stopwatch_pause(t_sw_cr)

    !$acc host_data use_device(zmel, zsec)
    ierr = gemm(czmelwc, zmel, zsec, ntqxx, ntqxx, nbb*(ns2-ns1+1), opA = m_op_T, beta = CONE, ldC = ntq)
    !$acc end host_data
    !$acc kernels loop independent
    do itp = 1, ntqxx
      zsec(itp,itp) = dreal(zsec(itp,itp))+img*min(dimag(zsec(itp,itp)),0d0) !enforce Imzsec<0
    enddo
    !$acc end kernels
    deallocate(wv, wc, czmelwc)
  end subroutine get_correlation
  subroutine setwv()
    integer :: iqini, iqend, iw
    character(10) :: i2char
    real(8), parameter :: gb = 1000*1000*1000
    complex(kind=kp), allocatable :: wv(:,:)
    if(allocated(wvi)) deallocate(wvi)
    if(allocated(wvr)) deallocate(wvr)

    if(.not.any(kx==kxc(:))) return
    open(newunit=ifrcwi,file='WVI.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
    open(newunit=ifrcw, file='WVR.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)

    if(.not.keepwv) return
    write(stdo, ftox) 'save WVI and WVR on CPU or GPU (if GPU is used) memory. This requires sufficient memory'
    call flush(stdo)
    allocate(wv(nblochpmx, nblochpmx), wvi(nblochpmx,nblochpmx,niw))
    do iw = 1, niw
      read(ifrcwi,rec=iw) wv(:,:)
      wvi(:,:,iw) = wv !send to GPU
    enddo
    allocate(wvr(nblochpmx,nblochpmx,nw_i:nw))
    do iw = nw_i, nw
      read(ifrcw,rec=iw-nw_i+1) wv
      wvr(:,:,iw) = wv
    enddo
    write(stdo, '(X,A,2F8.3)') 'WVI/WVR : sizes (GB)', dble(size(wvi))*kp*2/gb, dble(size(wvr))*kp*2/gb
    deallocate(wv)
  end subroutine setwv
  subroutine releasewv()
    if(.not.any(kx==kxc(:))) return
    if(allocated(wvi)) deallocate(wvi)
    if(allocated(wvr)) deallocate(wvr)
    close(ifrcwi)
    close(ifrcw)
  end subroutine releasewv
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
endmodule m_sxcf_gpu
