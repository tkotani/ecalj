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
!!    eq      = eigenvalues at q
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
!!    nlmto   = total number of MTO+lo basis functions
!!    ngrp    = no. group elements (rotation matrices)
!!    niw     = no. frequencies along the imaginary axis
!!    nw_i:nw  = no. frequencies along the real axis. nw_i=0 or -nw.
!!    zsec(itp,itpp,iq)> = <psi(itp,q(:,iq)) |SEc| psi(iq,q(:,iq)>
!! \endverbatim
module m_sxcf_main
  use m_readqg,only   : Readqg0
  use m_readeigen,only: Readeval
  use m_keyvalue,only   : Getkeyvalue
  use m_zmel,only : Get_zmel_init, get_zmel_init_gpu, Setppovlz,  zmel !Setppovlz_chipm,
  use m_itq,only: itq,ntq,nbandmx
  use m_genallcf_v3,only: nlmto,nspin,nctot,niw,ecore !,symgg
  use m_read_bzdata,only: qibz,qbz,wk=>wbz,nqibz,nqbz,wklm,lxklm,nq0i, wqt=>wt,q0i, irk
  use m_readVcoud,only:   Readvcoud, vcoud,ngb,ngc
  use m_readfreq_r,only: freq_r, nw_i,nw,freqx,wx=>wwx,nblochpmx,mrecl,expa_,npm,nprecx
  use m_rdpp,only: Rdpp, nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
  use m_readqg,only:  ngpmx,ngcmx
  use m_readhbe,only: nband,mrecg
  use m_readgwinput,only: ua_,  corehole,wcorehole
  use m_ftox
  use m_lgunit,only:stdo
  use m_sxcf_count,only: ncount,kxc,nstateMax,nstti,nstte,nstte2, nwxic, nwxc,icountini,icountend,irkip
  use m_nvfortran,only:findloc
  use m_hamindex,only: ngrp
  use m_data_gpu, only: SetDataGPU_inkx, ExitDataGPU_inkx, SetDataGPU, ExitDataGPU
  implicit none
  public sxcf_scz_main, zsecall
  complex(8),allocatable,target:: zsecall(:,:,:,:) !output
  private
contains
  subroutine sxcf_scz_main(ef,esmr,exchange,ixc,nspinmx) 
    intent(in)             ef,esmr,exchange,ixc,nspinmx
    logical :: exchange
    integer :: ip,it,itp,i,ix,kx,irot,kr,nt0p,nt0m,nstate,nbmax,ntqxx,nt,iw,ivc,ifvcoud,ngb0,ifwd,nrot,nwp,ierr,iqini,iqend
    integer :: invr,ia,nn,ntp0,no,itpp,nrec,itini,itend,nbmxe,iwp,nwxi,nwx,iir, igb1,igb2,ix0,isp,nspinmx
    integer :: invrot,nocc,nlmtobnd,nt0,verbose,ififr, istate, nt_max ,noccx,ns2r, kxold,nccc,icount0
    integer:: ixx,ixc,icount,ns1,ns2
    real(8) :: ekc(nctot+nband),ekq(nband),det,q(3),wtt,wfac,qvv(3),eq(nband),omega(ntq),quu(3),ef,esmr,&
         freqw,ratio,qibz_k(3),qbz_kr(3),vc,omega0,omg, polinta, wexx,ua2_(niw),freqw1,q_r(3),qk(3)
    logical:: tote=.false., iprx,cmdopt0, oncew, onceww !, eibz4sig  
    character(10) :: i2char
    real(8),parameter :: wfaccut=1d-8,tolq=1d-8, pi=4d0*datan(1d0), fpi=4d0*pi, tpi=8d0*datan(1d0),ddw=10d0
    complex(8), parameter :: img=(0d0,1d0)
    logical,parameter :: debug=.false.,timemix=.false.
    complex(8),allocatable,target:: zwz(:,:,:),zw(:,:)
    real(8),allocatable:: we_(:,:),wfac_(:,:)
    complex(8),allocatable:: w3p(:),wtff(:)
    complex(8),pointer::zsec(:,:), ww(:,:)
    integer,allocatable :: ifrcw(:),ifrcwi(:),ndiv(:),nstatei(:,:),nstatee(:,:)!,irkip(:,:,:,:)
    logical,parameter:: cache=.true.
    logical:: emptyrun, use_gpu
    emptyrun=cmdopt0('--emptyrun')
    use_gpu = cmdopt0('--gpu')
    if(nw_i/=0) call rx('Current version we assume nw_i=0. Time-reversal symmetry')
    iqini = 1
    iqend = nqibz             
    allocate(zsecall(ntq,ntq,nqibz,nspinmx),source=(0d0,0d0)) 
    if(.not.exchange.and.(.not.emptyrun)) then! Read WV* containing W-v in MPB
       allocate(ifrcw(iqini:iqend),ifrcwi(iqini:iqend))
       do kx=iqini,iqend
          if(any(kx==kxc(:))) then !only for requied files for this rank
             open(newunit=ifrcw(kx), file='WVR.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
             open(newunit=ifrcwi(kx),file='WVI.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
          endif
       enddo
     endif
    ! NOTE: sum for G\timesW is controlloed by irkip, icountini:icountend
    kxold=-9999  ! To make MAINicountloop 3030 as parallel loop, set cache=.false.
    kxloop:                  do kx  =1,nqibz   ! kx is irreducible !kx is main axis where we calculate W(kx).
      qibz_k = qibz(:,kx)
      call Readvcoud(qibz_k,kx,NoVcou=.false.)  !Readin ngc,ngb,vcoud ! Coulomb matrix
      call Setppovlz(qibz_k,matz=.true.,npr=ngb)        !Set ppovlz overlap matrix used in Get_zmel_init in m_zmel
      irotloop:              do irot=1,ngrp    ! (kx,irot) determines qbz(:,kr), which is in FBZ. W(kx) is rotated to be W(g(kx))
          iploopexternal:    do ip=1,nqibz     !external index for q of \Sigma(q,isp)
            isploopexternal: do isp=1,nspinmx  !external index
              kr = irkip(isp,kx,irot,ip)
              if(kr==0) cycle
              q(1:3)= qibz(1:3,ip)
              qbz_kr= qbz (:,kr)   !rotated qbz vector.
              qk =  q - qbz_kr     !<M(qbz_kr) phi(q-qbz_kr)|phi(q)>
              eq = readeval(q,isp) !readin eigenvalue
              omega(1:ntq) = eq(1:ntq)  !1:ntq
              ekq = readeval(qk, isp) 
              ekc(1:nctot)= ecore(1:nctot,isp) ! core
              ekc(nctot+1:nctot+nband) = ekq (1:nband)
              nt0  = count(ekc<ef) 
              nt0p = count(ekq<ef+ddw*esmr) +nctot 
              nt0m = count(ekq<ef-ddw*esmr) +nctot
              wtt = wk(kr)
              ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
              zsec => zsecall(1:ntqxx,1:ntqxx,ip,isp)
            NMBATCHloop: do 3030 icount = icountini(isp,ip,irot,kx),icountend(isp,ip,irot,kx) !batch of middle states.
              call flush(stdo)
              ns1 =nstti(icount) ! Range of middle states is [ns1:ns2] for given icount
              ns2 =nstte(icount) ! 
              nwxi=nwxic(icount) !minimum omega for W
              nwx =nwxc(icount)  !max omega for W
              ns2r=nstte2(icount) !Range of middle states [ns1:ns2r] for CorrelationSelfEnergyRealAxis
              ZmelBlock:block !zmel(ib=ngb,it=ns1:ns2,itpp=1:ntqxx)= <M(qbz_kr,ib) phi(it,q-qbz_kr,isp) |phi(itpp,q,isp)> 
                if(emptyrun) goto 1212
                if(use_gpu) then
                  call get_zmel_init_gpu(q,qibz_k,irot,qbz_kr, ns1,ns2,isp, 1,ntqxx, isp,nctot,ncc=0,iprx=debug, &
                                         zmelconjg=.false.)! Get zmel(ngb,ns1:ns2,ntqxx)
                else
                  call get_zmel_init(q,qibz_k,irot,qbz_kr, ns1,ns2,isp, 1,ntqxx, isp,nctot,ncc=0,iprx=debug,zmelconjg=.false.)! Get zmel(ngb,ns1:ns2,ntqxx)
                endif
1212            continue 
              endblock ZmelBlock
              ExchangeMode: if(exchange) then      
                if(use_gpu) then
                  ExchangeSelfEnergyGPU: Block
                    use m_xc_gpu, only: get_exchange
                    call get_exchange(kx, isp, ef, ekc, esmr, ns1, ns2, ntqxx, wtt, zsecall(1,1,ip,isp), &
                                      emptyrun = emptyrun)
                  EndBlock ExchangeSelfEnergyGPU
                else
                ExchangeSelfEnergy: Block
                  real(8):: wfacx,vcoud_(ngb),wtff(ns1:ns2) !range of middle states ns1:ns2
                  !character(8):: xt ;call timeshow("ExchangeMODE1 icount="//trim(xt(icount)))
                  vcoud_= vcoud                                    ! kx==1 must be for q=0     
                  if(kx==1) vcoud_(1)=wklm(1)*fpi*sqrt(fpi)/wk(kx) ! voud_(1) is effective v(q=0) in the Gamma cell. 
                  ! wtff = [(1d0,it=ns1,nctot), (wfacx(-1d99, ef, ekc(it), esmr),it=max(nctot+1,ns1),ns2)] !bugfix 2023-5-18 ns1==>max(nctot+1,ns1)
                  wtff(ns1:nctot) =1d0 !these are for nvfortran24.1
                  do it=max(nctot+1,ns1),ns2
                    wtff(it) = wfacx(-1d99, ef, ekc(it), esmr) !bugfix 2023-5-18 ns1==>max(nctot+1,ns1)
                  enddo
                  if(corehole) wtff(ns1:nctot) = wtff(ns1:nctot) * wcorehole(ns1:nctot,isp)
                  do concurrent(itp=1:ntqxx, itpp=1:ntqxx)
                    if(emptyrun) cycle !probably not so slow but for no error for --emptyrun
                    zsec(itp,itpp)=zsec(itp,itpp) - wtt* &
                         sum( [(sum(dconjg(zmel(:,it,itp))*vcoud_(:)*zmel(:,it,itpp))*wtff(it),it=ns1,ns2)] ) !this may work even for nvfortran24.1
                  enddo
                EndBlock ExchangeSelfEnergy
                endif
                if(timemix) call timeshow("ExchangeSelfEnergy cycle")
                cycle  
              endif ExchangeMode
              if(use_gpu) then
                !$acc update host(zmel)
              endif
              CorrelationMode: Block! See Eq.(55) around of PRB76,165106 (2007) !range of middle states is [ns1:ns2]
                integer:: nm
                complex(8):: zmelc  (1:ntqxx,ns1:ns2,1:ngb)
                complex(8):: zmelcww(1:ntqxx,ns1:ns2,1:ngb)
                if(.not.emptyrun) zmelc = reshape(dconjg(zmel),shape=shape(zmelc),order=[3,2,1]) !notslow. but for no error for --emptyrun 
                if(timemix) call timeshow(" CorrelationSelfEnergyImagAxis:")
                CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)! Integration along ImAxis for zwz(omega) 
                  use m_readfreq_r,only: wt=>wwx,x=>freqx
                  real(8),parameter:: rmax=2d0
                  real(8):: wgtim_(0:npm*niw),wgtim(0:npm*niw,ntqxx,ns1:ns2),we,cons(niw),omd(niw),omd2w(niw)
                  complex(8),target:: zw(nblochpmx,nblochpmx)!,zwz(ns1:ns2,ntqxx,ntqxx)
                  itpitdo:do concurrent(itp =1:ntqxx, it=ns1:ns2)
                    we=.5d0*(omega(itp)-ekc(it)) !we in hartree unit (atomic unit)
                    associate(sig=>.5d0*esmr, sig2=>2d0*(.5d0*esmr)**2, aw=>abs(ua_*we), aw2=>(ua_*we)**2)
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
                    endassociate
                    wgtim(:,itp,it)=wtt*wgtim_ !! Integration weight wgtim along im axis for zwz(0:niw*npm)
                  enddo itpitdo
                  zmelcww=0d0
                  iwimag:do ixx=0,niw !niw is ~10. ixx=0 is for omega=0 nw_i=0 (Time reversal) or nw_i =-nw
                    if(emptyrun) cycle
                    !   iwimag:do concurrent(ixx=0:niw)  !concurrent may(or maynot) be problematic because zmelcww is for reduction 
                    if(ixx==0)read(ifrcw(kx),rec=1+(0-nw_i)) zw !direct access Wc(0) = W(0)-v ! nw_i=0 (Time reversal) or nw_i =-nw
                    if(ixx>0) read(ifrcwi(kx),rec=ixx) zw       ! direct access read Wc(i*omega)=W(i*omega)-v
                    ww=>zw(1:ngb,1:ngb)
                    call matmaw(zmelc,ww,zmelcww, ntqxx*(ns2-ns1+1),ngb,ngb, wgtim(ixx,:,:)) ! time-consuming
                  enddo iwimag
                EndBlock CorrelationSelfEnergyImagAxis
                CorrelationSelfEnergyRealAxis: Block !Real Axis integral. Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
                  use m_wfac,only: wfacx2, weavx2
                  integer:: itini,itend
                  integer:: iwgt3(ns1:ns2r,ntqxx),i1,i2,iw,ikeep,ix,ixs,nit_(ntqxx,nwxi:nwx),icountp,ncoumx,iit,irs
                  real(8):: we_(ns1:ns2r,ntqxx),wfac_(ns1:ns2r,ntqxx),omg,esmrxx,wgt3(0:2,ns1:ns2r,ntqxx),amat(3,3)
                  !3-point interpolation weight for we_(it,itp) 
                  complex(8),target:: zw(nblochpmx,nblochpmx)
                  complex(8):: zadd(ntqxx),wv3(ngb,ngb,0:2)
                  logical:: itc(ns1:ns2r,nwxi:nwx,ntqxx)
                  if(timemix) call timeshow(" CorrelationSelfEnergyRealAxis:")
                  itc=.false.
                  do concurrent(itp=1:ntqxx)
                    omg  = omega(itp)
                    itini= merge(max(ns1,nt0m+1),  ns1, mask= omg>=ef)
                    itend= merge(ns2r,  min(nt0p,ns2r), mask= omg>=ef)
                    do concurrent(it=itini:itend)     ! nt0p corresponds to efp
                      wfac_(it,itp) = wfacx2(omg,ef, ekc(it), merge(0d0,esmr,mask=it<=nctot)) !Gaussian smearing 
                      if(wfac_(it,itp)<wfaccut) cycle 
                      wfac_(it,itp)=  wfac_(it,itp)*wtt*dsign(1d0,omg-ef) !wfac_ = $w$ weight (smeared thus truncated by ef). See the sentences.
                      we_(it,itp)  = .5d0*abs(omg-weavx2(omg,ef, ekc(it),esmr)) !we_= \bar{\omega_\epsilon} in sentences next to Eq.58 in PRB76,165106 (2007)
                      ixs = findloc(freq_r(1:nw)>we_(it,itp),value=.true.,dim=1)
                      iwgt3(it,itp) = ixs+1-2    !Starting omega index ix for it,itp    ! iwgt3(it,itp) = iirx(itp)*(ixs+1-2) !iirx(ntqxx),
                      itc(it,ixs+1-2,itp)=.true. !W(we_(it,itp))=\sum_{i=0}^2 W(:,:,ix+i)*wgt3(i)         
                      associate(x=>we_(it,itp),xi=>freq_r(ixs-1:ixs+1))!x=>we_ is \omega_\epsilon in Eq.(55). 
                        amat(1:3,1)= 1d0                   !old version: call alagr3z2wgt(we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
                        amat(1:3,2)= xi(1:3)**2
                        amat(1:3,3)= xi(1:3)**4
                        wgt3(:,it,itp)= wfac_(it,itp)*matmul([1d0,x**2,x**4], inverse33(amat)) 
                      endassociate
                    enddo
                  enddo
                  !iirx=1 !where(omega<ef.and.nw_i/=0) iirx=-1 !if(any(iirx(1:ntqxx)/=1))call rx('sxcf: iirx=-1(TR breaking) is not yet implemented')
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
              EndBlock CorrelationMode
              if(timemix) call timeshow("   end icount do 3030")
3030        enddo NMBATCHloop
          enddo isploopexternal
        enddo iploopexternal
      enddo irotloop
      call ExitDataGPU_inkx()
    enddo kxloop
    call ExitDataGPU()
  write(stdo,ftox)'endof 3030loop'
  endsubroutine sxcf_scz_main
  pure subroutine matmaw(a,b,c,n1,n2,n3,ww)
    integer, intent(in) :: n1,n2,n3
    complex(8), intent(in) :: a(n1,n2), b(n2,n3)
    real(8), intent(in) :: ww(n1)
    complex(8), intent(inout) :: c(n1,n3)
    complex(8):: aa(n1,n2)
    integer:: ix
    do ix=1,n2
       aa(:,ix) = ww(:)*a(:,ix)
    enddo
    c= c+ matmul(aa,b)
  end subroutine matmaw
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
endmodule m_sxcf_main
