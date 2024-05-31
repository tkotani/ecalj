module m_xc_gpu
  use m_blas, only: m_op_c, m_op_n, m_op_t, zmm, zmm_batch
  use m_zmel, only: zmel, nbb  ! size (nbb,ns1:ns2,nqtot)
  use m_readVcoud,only: vcoud, ngb
  use m_genallcf_v3, only: nctot, niw
  use m_itq, only: ntq
  use m_readfreq_r, only: freq_r, nw_i, nw, freqx, wx=>wwx, nblochpmx, expa_, npm
  use m_readhbe, only: nband
  use m_readgwinput, only: ua_, corehole, wcorehole
  logical, parameter:: cache=.true.
  complex(8), allocatable :: wc_k(:,:,:)
  integer :: k_of_wc = 0
  logical :: has_wc_k = .false.
contains
  ! subroutine set_wc_k(kx, ifrcwkx, ifrcwikx)
  !   integer, intent(in) :: kx, ifrcwkx, ifrcwikx
  !   if(allocated(wc_k)) deallocate(wc_k)
  !   allocate(wc_k(
  !     iwimag:do ixx = 1, nw
  !       if(ixx==0)read(ifrcw(kx),rec=1+(0-nw_i)) zw(:,:,iw) !direct access Wc(0) = W(0)-v ! nw_i=0 (Time reversal) or nw_i =-nw
  !       if(ixx>0) read(ifrcwi(kx),rec=ixx) zw       ! direct access read Wc(i*omega)=W(i*omega)-v
  !     enddo iwimag
  ! end subroutine
  subroutine get_exchange(kx, isp, ef, ekc, esmr, ns1, ns2, ntqxx, wtt, zsec, emptyrun)
    use m_read_bzdata, only: wk=>wbz, wklm
    implicit none
    integer, intent(in) :: kx, isp, ns1, ns2, ntqxx
    real(8), intent(in) :: ef, esmr, wtt, ekc(:)
    complex(8), intent(inout) ::zsec(ntq,ntq)
    logical, intent(in), optional :: emptyrun
    real(8), parameter :: pi = 4d0*datan(1d0), fpi = 4d0*pi
    real(8) :: wfacx ! external function
    real(8), allocatable :: vcoud_buf(:), wtff(:)
    complex(8), allocatable :: vzw(:,:,:)
    integer :: it, itp, itpp, igb, ierr
#ifdef __GPU
    attributes(device) :: vzw, vcoud_buf
#endif
    if(present(emptyrun)) then
      if(emptyrun) return
    endif
    allocate(wtff(ns1:ns2))
    wtff(ns1:nctot) = 1d0
    do it = max(nctot+1,ns1), ns2
      wtff(it) = wfacx(-1d99, ef, ekc(it), esmr)
    enddo
    if(corehole) wtff(ns1:nctot) = wtff(ns1:nctot) * wcorehole(ns1:nctot,isp)

    allocate(vcoud_buf(ngb))
    !$acc data copyin(vcoud, wklm(1), wk(1), wtff) copy(zsec)
    !$acc kernels
    vcoud_buf(1:ngb) = vcoud(1:ngb)
    !$acc end kernels
    if(kx == 1) vcoud_buf(1) = wklm(1)*fpi*sqrt(fpi)/wk(1)
    allocate(vzw(ngb,ns1:ns2,ntqxx))
    !$acc kernels loop independent collapse(3)
    do itpp = 1, ntqxx
      do it = ns1, ns2 
        do igb = 1, ngb
          vzw(igb,it,itpp) = vcoud_buf(igb)*zmel(igb,it,itpp)*wtff(it)
        enddo
      enddo
    enddo
    !$acc end kernels
    !$acc host_data use_device(zmel)
    ierr = zmm(zmel(1,ns1,1), vzw, zsec, ntqxx, ntqxx, (ns2-ns1+1)*ngb, opA = m_op_c, &
               alpha = dcmplx(-wtt,0D0), beta = (1D0,0d0), ldA = (ns2-ns1+1)*nbb, ldC = ntq)
    !$acc end host_data
    !$acc end data
    deallocate(vzw, vcoud_buf, wtff)
  end subroutine
#ifdef __SKIP
  subroutine get_correlation(ns1, ns2, ntqxx, nt0p, nt0m, ifrcwkx, ifrcwikx, ekc, omega, wtt, emptyrun)
    integer, intent(in) :: ns1, ns2, ntqxx, nt0p, nt0m, ifrcwkx, ifrcwikx
    real(8), intent(in) :: omega(ntq), ekc(:), wtt
    complex(8), allocatable :: zmelc(:,:,:), zmelcww(:,:,:)
    logical, intent(in), optional :: emptyrun
    complex(8), allocatable :: zw(:,:)
    #ifdef __GPU
      attributes(device) :: zmelw, zmelwzw
    #endif
    integer:: nm
     ! zmel(nbb,ns1:ns2,nqtot)
     ! zmelc(1:ntqxx,ns1:ns2,1:ngb)
    allocate(zw(nblochpmx,nblochpmx))
    allocate(zmelw(1:ngb,ns1:ns2,1:ntqxx))
    allocate(zmelwzw(1:ngb,ns1:ns2,1:ntqxx))
    if(timemix) call timeshow(" CorrelationSelfEnergyImagAxis:")
    CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)! Integration along ImAxis for zwz(omega) 
      use m_readfreq_r,only: wt=>wwx,x=>freqx
      real(8),parameter:: rmax=2d0
      real(8):: wgtim_(0:npm*niw), wgtim(0:npm*niw,ntqxx,ns1:ns2), we, cons(niw), omd(niw), omd2w(niw)
      real(8):: sig = .5d0*esmr, sig2=2d0*(.5d0*esmr)**2
      real(8):: aw, aw2
      !$acc data create(wgtim_, wgtim, omd, omd2w) copyin(omega(1:ntqxx),ekc(ns1:ns2))
      !$acc kernels loop independent collapse(2)
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
          wgtim(:,itp,it)=wtt*wgtim_ !! Integration weight wgtim along im axis for zwz(0:niw*npm)
        enddo 
      enddo 
      !$acc end kernels
      zmelwzw = (0d0, 0D0)
      iwimag:do ixx = 0, niw !niw is ~10. ixx=0 is for omega=0 nw_i=0 (Time reversal) or nw_i =-nw
        if(emptyrun) cycle
        if(ixx==0)read(ifrcw(kx),rec=1+(0-nw_i)) zw !direct access Wc(0) = W(0)-v ! nw_i=0 (Time reversal) or nw_i =-nw
        if(ixx>0) read(ifrcwi(kx),rec=ixx) zw       ! direct access read Wc(i*omega)=W(i*omega)-v
        !$acc kernels loop independent collapse(2)
        do itp = 1, ntqxx
          do it = ns1, ns2
            zmelw(1:ngb,it,itp) = zmel(1:ngb,it,itp)*wgtim(ixx,itp,it) 
          enddo
        enddo
        !$acc end kernels
        !$acc data copyin(zw)
        ierr = zmm(zmelw, zw, zmelwzw, ngb, ntqxx*(ns2-ns1+1), ngb, opA = m_op_c, beta = (1D0, 0D0), ldB = nblochpmx)   
        !$acc end data
      enddo iwimag
    EndBlock CorrelationSelfEnergyImagAxis
    CorrelationSelfEnergyRealAxis: Block !Real Axis integral. Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
      use m_wfac, only: wfacx2, weavx2
      integer :: itini,itend
      integer :: iwgt3(ns1:ns2r,ntqxx),iw,ikeep,ix,ixs,nit_(ntqxx,nwxi:nwx),irs
      real(8) :: we_(ns1:ns2r,ntqxx),wfac_(ns1:ns2r,ntqxx),omg,esmrxx,wgt3(0:2,ns1:ns2r,ntqxx),amat(3,3)
      !3-point interpolation weight for we_(it,itp) 
      complex(8):: zadd(ntqxx),wv3(ngb,ngb,0:2)
      logical:: itc(ns1:ns2r,nwxi:nwx,ntqxx)
      if(timemix) call timeshow(" CorrelationSelfEnergyRealAxis:")
      itc=.false.

      !$acc kernels
      do itp=1, ntqxx
        omg  = omega(itp)
        itini= merge(max(ns1,nt0m+1),  ns1, mask= omg>=ef)
        itend= merge(ns2r,  min(nt0p,ns2r), mask= omg>=ef)
        do it=itini, itend     ! nt0p corresponds to efp
          wfac_(it,itp) = wfacx2(omg,ef, ekc(it), merge(0d0,esmr,mask=it<=nctot)) !Gaussian smearing 
          if(wfac_(it,itp)<wfaccut) cycle 
          wfac_(it,itp)=  wfac_(it,itp)*wtt*dsign(1d0,omg-ef) !wfac_ = $w$ weight (smeared thus truncated by ef). See the sentences.
          we_(it,itp)  = .5d0*abs(omg-weavx2(omg,ef, ekc(it),esmr)) !we_= \bar{\omega_\epsilon} in sentences next to Eq.58 in PRB76,165106 (2007)
          ixs = findloc(freq_r(1:nw)>we_(it,itp),value=.true.,dim=1)
          iwgt3(it,itp) = ixs+1-2    !Starting omega index ix for it,itp    ! iwgt3(it,itp) = iirx(itp)*(ixs+1-2) !iirx(ntqxx),
          itc(it,ixs+1-2,itp)=.true. !W(we_(it,itp))=\sum_{i=0}^2 W(:,:,ix+i)*wgt3(i)         
          x = we_(it,itp)
          xi(1:3)=freq_r(ixs-1:ixs+1) !x=>we_ is \omega_\epsilon in Eq.(55). 
          amat(1:3,1)= 1d0                   !old version: call alagr3z2wgt(we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
          amat(1:3,2)= xi(1:3)**2
          amat(1:3,3)= xi(1:3)**4
          wgt3(:,it,itp)= wfac_(it,itp)*matmul([1d0,x**2,x**4], inverse33(amat)) 
        enddo
      enddo
      !$acc end kernels

      ikeep=99999
      do ix = nwxi,nwx  
        if(all(.not.itc(:,ix,:))) cycle 
        irs=0
        if(cache) then ! use wv3 at previous ix. a cache mechanism, notefficient?
          if(ikeep+1==ix) then 
            wv3(:,:,0)=wv3(:,:,1) 
            wv3(:,:,1)=wv3(:,:,2)
            irs=2
          elseif(ikeep+2==ix) then 
            wv3(:,:,0)=wv3(:,:,2) 
            irs=1
          endif
        endif
        if(emptyrun) cycle
        do iw=irs,2 !Set wv3(:,:,0:2). 0:2 means ix:ix+2
          read(ifrcw(kx),rec=iw+ix-nw_i+1) zw ! direct access Wc(omega) = W(omega) - v
          ww=>zw(1:ngb,1:ngb)
          wv3(:,:,iw)=(ww+transpose(dconjg(ww)))/2d0 !hermite part
        enddo ! wv3 should contain zw(ix+0:ix+2) now
        ikeep=ix
        do itp=1,ntqxx !lbound(zsec,1),ubound(zsec,1)
          do it =ns1,ns2r ! for it for given itp,ix
            if(.not.itc(it,ix,itp)) cycle  !wv33 gives interpolated value of W(we_(it,itp))
            zmelcww(itp,it,:) = zmelcww(itp,it,:) +  matmul(zmelc(itp,it,:), &
                 wv3(:,:,0)*wgt3(0,it,itp) +wv3(:,:,1)*wgt3(1,it,itp) +wv3(:,:,2)*wgt3(2,it,itp) ) !time-consuming
          enddo
        enddo

      enddo

      if(timemix) call timeshow(" End of CorrelationSelfEnergyRealAxis:")
    EndBlock CorrelationSelfEnergyRealAxis
    nm = (ns2-ns1+1)*ngb
    if(emptyrun) return
    zsec= zsec+ matmul(reshape(zmelcww,[ntqxx,nm]),reshape(&
         reshape(zmel,shape=[ns2-ns1+1,ngb,ntqxx],order=[2,1,3]), [nm,ntqxx])) !time-consuming probably.
    forall(itp=1:ntqxx) zsec(itp,itp)=dreal(zsec(itp,itp))+img*min(dimag(zsec(itp,itp)),0d0) !enforce Imzsec<0
  end subroutine
#endif
end module
