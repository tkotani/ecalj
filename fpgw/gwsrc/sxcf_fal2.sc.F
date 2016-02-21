!> this module is only because name=name argument binding. No data
      module m_sxcfsc
      contains
      subroutine sxcf_fal3_scz(kount,qip,itq,ntq,ef,esmr,
     i nsp,isp,             
     i qbas,ginv,
     i qibz,qbz,wk,nstbz,irkip,nrkip, 
     i freq_r,nw_i,nw, freqx,wx,dw,         
     i ecore,
     i nlmto,nqibz,nqbz,nctot,
     i nbloch,ngrp,niw,nq,     
     i nblochpmx ,ngpmx,ngcmx, 
     i wgt0,nq0i,q0i,symgg, alat, nband, ifvcfpout,
     i exchange,screen,cohtest, ifexsp,
     i nbmx,ebmx,           
     i wklm,lxklm,          
     i eftrue,
     i jobsw,                   != iSigma_en 
     i hermitianW,
     o zsec,coh,nbandmx)
      use m_readqg,only   : readqg0
      use m_readeigen,only: readeval
      use m_keyvalue,only   : getkeyvalue
      use m_zmel,only     : get_zmelt,
     i ppovlz, 
     o zmel,zmeltt 
      implicit none
!> \brief
!! Calcualte full simga_ij(e_i)= <i|Re[Sigma](e_i)|j> 
!! ---------------------
!! \param exchange 
!!   - T : Calculate the exchange self-energy
!!   - F : Calculate correlated part of the self-energy
!! \param zsec
!!   - S_ij= <i|Re[S](e_i)|j>
!!   - Note that S_ij itself is not Hermite becasue it includes e_i.
!!     i and j are band indexes
!! \param coh dummy 
!! \param screen dummy 
!!
!! \remark
!! \verbatim
!! Jan2013: eftrue is added.
!!   ef=eftrue(true fermi energy) for valence exchange and correlation mode.
!!   but ef is not the true fermi energy for core-exchange mode.
!!
!! Jan2006
!!     "zsec from im-axis integral part"  had been symmetrized as
!!     &        wtt*.5d0*(   sum(zwzi(:,itp,itpp))+ !S_{ij}(e_i)
!!     &        dconjg( sum(zwzi(:,itpp,itp)) )   ) !S_{ji}^*(e_j)= S_{ij}(e_j)
!!     However, I now do it just the 1st term.
!!     &        wtt* sum(zwzi(:,itp,itpp))   !S_{ij}(e_i)
!!     This is OK because the symmetrization is in hqpe.sc.F
!!     Now zsec given in this routine is simply written as <i|Re[S](e_i)|j>.
!!     ( In the version until Jan2006 (fpgw032f8), only the im-axis part was symmetrized.
!!     But it was not necessary from the begining because it was done in hqpe.sc.F
!!     
!!     (Be careful as for the difference between
!!     <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
!!     ---because e_i is included.
!!     The symmetrization (hermitian) procedure is inlucded in hqpe.sc.F
!!
!!     NOTE: matrix element is given by "call get_zmelt". It returns  zmelt or zmeltt.
!!
!! jobsw switch
!!  1-5 scGW mode.
!!   diag+@EF      jobsw==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
!!   xxx modeB (Not Available now)  jobsw==2 SE_nn'((e_n+e_n')/2)  !we need to recover comment out for jobsw==2, and test.
!!   mode A        jobsw==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2 (Usually usued in QSGW).
!!   @Ef           jobsw==4 SE_nn'(ef) 
!!   diagonly      jobsw==5 delta_nn' SE_nn(e_n) (not efficient memoryuse; but we don't use this mode so often).
!!
!! Output file in hsfp0 should contain hermitean part of SE
!!    ( hermitean of SE_nn'(e_n) means SE_n'n(e_n')^* )
!!             we use that zwz(itp,itpp)=dconjg( zwz(itpp,itp) )
!! Caution! npm=2 is not examined enough...
!!
!! Calculate the exchange part and the correlated part of self-energy.
!! T.Kotani started development after the analysis of F.Aryasetiawan's LMTO-ASA-GW.
!! We still use some of his ideas in this code.
!!
!! See paper   
!! [1]T. Kotani and M. van Schilfgaarde, ??Quasiparticle self-consistent GW method: 
!!     A basis for the independent-particle approximation, Phys. Rev. B, vol. 76, no. 16, p. 165106[24pages], Oct. 2007.
!! [2]T. Kotani, Quasiparticle Self-Consistent GW Method Based on the Augmented Plane-Wave 
!!    and Muffin-Tin Orbital Method, J. Phys. Soc. Jpn., vol. 83, no. 9, p. 094711 [11 Pages], Sep. 2014.
!!
!! =======================================================================
!! Omega integral for SEc
!!   The integral path is deformed along the imaginary-axis, but together with contribution of poles.
!!   See Fig.1 and around in Ref.[1].
!!
!! ---Integration along imaginary axis.---
!!   ( Current version for it, wintzsg_npm, do not assume time-reversal when npm=2.)
!!   Integration along the imaginary axis: -----------------
!!    (Here is a memo by F.Aryasetiawan.)
!!     (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
!!    Gaussian integral along the imaginary axis.  
!!    transform: x = 1/(1+w')
!!     this leads to a denser mesh in w' around 0 for equal mesh x
!!    which is desirable since Wc and the lorentzian are peaked around w'=0
!!     wint = - (1/pi) < [x=0,1] Wc(iw') (w-e)x^2/{(w-e)^2 + w'^2} >
!!     
!!     the integrand is peaked around w'=0 or x=1 when w=e
!!     to handel the problem, add and substract the singular part as follows:
!!     wint = - (1/pi) < [x=0,1] { Wc(iw') - Wc(0)exp(-a^2 w'^2) }
!!     * (w-e)/{(w-e)^2 +w'^2}x^2 >
!!     - (1/2) Wc(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
!!     
!!     the second term of the integral can be done analytically, which
!!     results in the last term a is some constant
!!     
!!     when w = e, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w') and
!!     the integral becomes -Wc(0)/2
!!     this together with the contribution from the pole of G (s.u.)
!!     gives the so called static screened exchange -Wc(0)
!!
!! ---Integration along real axis (contribution from the poles of G: SEc(pole))
!!    See Eq.(34),(55), and (58) and around in Ref.[1]. We now use Gaussian Smearing.
!! =======================================================================
!
!! ----------------------------------------------
!!     q     =qip(:,iq)  = q-vector in SEc(q,t). 
!!    itq     = states t at q
!!    ntq     = no. states t
!!    eq      = eigenvalues at q
!!     ef      = fermi level in Rydberg
!!   WVI, WVR: direct access files for W. along im axis (WVI) or along real axis (WVR)
!!   freq_r(nw_i:nw)   = frequencies along real axis. freq_r(0)=0d0
!!
!!    qbas    = base reciprocal lattice vectors
!!    ginv    = inverse of qbas s. indxrk.f
!!
!!     wk     = weight for each k-point in the FBZ
!!    qbz     = k-points in the 1st BZ
!!
!!    wx      = weights at gaussian points x between (0,1)
!!     ua_      = constant in exp(-ua^2 w'^2) s. wint.f
!!     expa    = exp(-ua^2 w'^2) s. wint.f
!!
!!    irkip(k,R,nq) = gives index in the FBZ with k{IBZ, R=rotation
!!
!!   nqibz   = number of k-points in the irreducible BZ
!!   nqbz    =                           full BZ
!!    natom   = number of atoms
!!    nctot   = total no. of allowed core states
!!    nbloch  = total number of Bloch basis functions
!!    nlmto   = total number of MTO+lo basis functions
!!    ngrp    = no. group elements (rotation matrices)
!!    niw     = no. frequencies along the imaginary axis
!!    nw_i:nw  = no. frequencies along the real axis. nw_i=0 or -nw.
!!    zsec(itp,itpp,iq)> = <psi(itp,q(:,iq)) |SEc| psi(iq,q(:,iq)>
!!
!! ----------------------------------------------
!! \endverbatim
      integer:: dummy4doxygen

! input variables
      logical, intent(in) :: exchange,screen,cohtest
      integer, intent(in) :: ntq,nqbz,nqibz,ngrp,nq,niw !,natom
      integer, intent(in) :: nband,nlmto,nq0i,nctot,isp,nsp !,mdim(*) !,nlnmx
      integer, intent(in) :: ifvcfpout,nbloch,nblochpmx !nl,nnc, nclass
      integer, intent(in) :: itq(ntq) !,nstar(nqibz) !miat(natom,ngrp),mdimx,
      integer, intent(in) :: irkip(nqibz,ngrp,nq),nrkip(nqibz,ngrp,nq)
      integer, intent(in) :: kount(nqibz,nq),ngpmx,ngcmx,ifexsp,jobsw
      integer, intent(in) :: nbmx(2) !,nlnmv(*),nlnmc(*)!,iclass(*),icore(*)
      integer, intent(in) :: nstbz(nqbz) !,nomx !,nkpo,nnmx,imdim(*)ncore(*),
      integer, intent(in) :: lxklm !,invg(ngrp) !nnr(nkpo),nor(nkpo),
c      integer, intent(in) :: il(*),in(*),im(*),nn_,lx(*),nx_(*),nxx_ !,nlnm(*)
      real(8), intent(in) :: wgt0(nq0i,ngrp),symgg(3,3,ngrp)
      real(8), intent(in) :: q0i(1:3,1:nq0i),alat,ecore(nctot) !shtvg(3,ngrp),
      real(8), intent(in) :: qbas(3,3),ginv(3,3)
      real(8), intent(in) :: wk(nqbz),qibz(3,nqibz) !tiat(3,natom,ngrp),
      real(8), intent(in) :: qbz(3,nqbz),freqx(niw),wx(niw),ef,esmr,dw
      real(8), intent(in) :: ebmx(2),wklm((lxklm+1)**2) !,qrr(3,nkpo)
      real(8), intent(in) :: qip(3,nq),eftrue

c      integer,intent(in):: iwini,iwend
c      real(8),optional::exx

! output variables
c      real(8),intent(in),optional:: freqsig(iwini:iwend)
      integer, intent(in) ::nbandmx(nq)
      complex(8), intent(out),optional :: zsec(ntq,ntq,nq) , coh(ntq,nq)
c      complex(8), intent(out),optional :: zsecd(iwini:iwend,ntq,nq)

! local variables
c     complex(8) :: zsecx(ntq,ntq,nq)          
c      complex(8), intent(in) :: pomatr(nnmx,nomx,nkpo)
c$$$      logical :: ua_auto !fixed to be .false.
c      real(8)::ppbrd ( 0:nl-1, nn_, 0:nl-1,nn_, 0:2*(nl-1),1:nxx_, 1:nsp*nclass)

      integer :: ifrcw,ifrcwi
      logical :: initp=.true.
      real(8),allocatable:: vcoud(:)

      integer :: ip, it, itp, i, ix, kx, irot, kr
      integer :: nt0p, nt0m,nstate , nbmax, ntqxx !iatomp(natom), 
      integer :: nt,nw,ixs,iw,ivc,ifvcoud,ngb0
      integer :: nprecx,mrecl,ifwd,nrot,nwp,nw_i,ierr 
      integer :: nstatetot,iqini,iqend, ngb,ngc !nbcut,
      integer :: invr,nbmxe,ia,nn,ntp0,no,itpp,nrec,npm,itini,itend
      integer :: iwp,nwxi,nwx,iir, igb1,igb2,ix0,iii

      real(8) :: tpi, ekc(nctot+nband),ekq(nband), det, q(3),ua_
      real(8) :: expa_(niw), qxx(3), symope(3,3),shtv(3) !tr(3,natom), 
      real(8) :: efp,efm,wtt,wfac,we,esmrx,qbasinv(3,3)
      real(8) :: qvv(3),pi,fpi,eq(nband),omega(ntq),quu(3),freqw,ratio
      real(8) :: qibz_k(3),qbz_kr(3),ddw,vc,omega0,omg

      complex(8) :: cphiq(nlmto,nband), cphikq(nlmto,nband)
      complex(8) :: zwzs0,zz2,zwz3(3)

! local arrays
      real(8),intent(in) :: freq_r(nw_i:nw)
      real(8),allocatable :: drealzzzmel(:,:,:), dimagzzzmel(:,:,:),uaa(:,:)
      complex(8),allocatable :: vcoul(:,:),w3p(:,:,:)   
      complex(8),allocatable :: zzzmel(:,:,:),zw (:,:)
      complex(8),allocatable :: zwz(:,:,:,:), zwz0(:,:,:),zwzi(:,:,:)
      complex(8),allocatable :: zwix(:,:),zwzix(:,:,:),zmel1(:) !,expikt(:)
      complex(8), allocatable :: zmel1_(:,:,:), zw3(:,:,:),zw3x(:,:)
      complex(8), allocatable :: zwz4(:,:),zwz44(:,:),pomat(:,:), zwzs(:)
      complex(8),allocatable :: ppovl(:,:),zcousq(:,:)
      complex(8),allocatable :: z1r(:,:),z2r(:,:),w3pi(:,:)

      real(8), parameter :: wfaccut=1d-8
      complex(8), parameter :: img=(0d0,1d0)

! external function
c      logical :: smbasis
c      logical :: test_symmetric_W
c      logical :: GaussSmear !fixed to be T
c      logical :: newaniso !fixed to be T
c      integer :: bzcase !fixed to be 1
      character(5) :: charnum5
      integer :: iopen,iclose
      integer :: invrot
      complex(8) ::  wintzsg_npm !wintzav,
      integer :: nocc
      real(8) :: wfacx
      real(8) :: wfacx2
      real(8) :: weavx2
      complex(8) :: alagr3z
      complex(8) :: alagr3z2

      integer:: ndummy1,ndummy2,nlmtobnd,nt0
      real(8):: wexx
c      complex(8),allocatable :: z1p(:,:,:),vcoult(:,:)
      logical :: debug, debugp,debug2=.false.
c      logical :: gass           !external
c      real(8):: wgtq0p
      integer::verbose,ififr,ifile_handle
      real(8):: ua2_(niw),freqw1
      integer :: istate,  nt_max !nbcutc,nbcutin, 
      real(8):: q_r(3),qk(3),omegat
      logical::  oncew, onceww, eibz4sig,  timemix

      integer,allocatable:: ixss(:,:),iirx(:)
      real(8),allocatable:: we_(:,:),wfac_(:,:)
      complex(8),allocatable:: zw3av(:,:),zmelw(:,:,:)
      integer:: noccx
      real(8)::polinta
      logical,allocatable:: ititpskip(:,:)

      logical:: tote=.false.
      logical:: hermitianW

      real(8),allocatable:: wcorehole(:,:)
      logical:: corehole
      integer:: ifcorehole
c      real(8),allocatable:: ppb(:)
c      allocate( ppb(nlnmx*nlnmx*mdimx*nclass))

c      real(8)::exxq

c--------------------------------------------------------------------
c!TIME0_0000 
c      write(6,*)'sxcf_fal3_scz'
      timemix=.false.
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      debug=.false.
      if(verbose()>=90) debug=.true.

cccccccccccccccc
      corehole=.false.
cccccccccccccccc

!! core-hole
      if(corehole) then
        ifcorehole=ifile_handle()
        open(ifcorehole,file='CoreHole') 
        if(allocated(wcorehole)) deallocate(wcorehole)
        allocate(wcorehole(nctot,nsp))
        do it=1,nctot
          read(ifcorehole,*) wcorehole(it,1:nsp)
        enddo
        close(ifcorehole)
      endif

      if(.not.exchange) then
        ifwd = iopen('WV.d',1,-1,0)
        read (ifwd,*) nprecx,mrecl
        ifwd = iclose('WV.d')
!! gauss_img : interpolation gaussion for W(i \omega).
        call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
        if(debug) write(6,*) ' sxcf_fal3_scz: Gausssmear=T'
        do ix = 1,niw           !! Energy mesh; along im axis.
          freqw     = (1d0 - freqx(ix))/ freqx(ix)
          expa_(ix) = exp(-(ua_*freqw)**2)
        enddo
        npm = 1                 ! npm=1    Timeveversal case
        if(nw_i/=0) npm = 2     ! npm=2 No TimeReversal case. Need negative energy part of W(omega)
      endif 

c      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
c      nbcutc=nctot+nbcut
      tpi         = 8d0*datan(1d0)
      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband


!!== ip loop to spedify external q ==
      do 1001 ip = 1,nq         
        if(sum(irkip(:,:,ip))==0) cycle ! next ip
        write (6,*) ip,'  out of ',nq,'  k-points(extrnal q) '
        q(1:3)= qip(1:3,ip)
        call readeval(q,isp,eq)
        do i  = 1,ntq
          omega(i) = eq(itq(i)) 
        enddo

!! we only consider bzcase()==1
        if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
        if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')

!! NOTE total number of 
!!    kx loop(do 1100) and irot loop (do 1000) makes all the k mesh points.
!!    When iqini=1 (Gamma point), we use effective W(q=0) defined in the paper.
        iqini=1
        iqend=nqibz             !no sum for offset-Gamma points.
        do 1100 kx = iqini,iqend 
          if(sum(irkip(kx,:,ip))==0) cycle ! next kx
!TIME0_01000
          write(6,*) ' ### do 1100 start kx=',kx,' from ',iqini,' through', iqend
c          if( kx <= nqibz ) then
            qibz_k= qibz(:,kx)
c          else
c            qibz_k= 0d0
c          endif
          if(timemix) call timeshow("11111 k-cycle")
          call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
          ngb = nbloch + ngc
          if(debug) write(6,*) ' sxcf: ngb=',ngb,nbloch

!! ===Readin diagonalized Coulomb interaction===
!!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
!!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
!!  Vould file is written in hvccfp0.m.F.
!! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
!! These can be also separeted into WVR.ID and WVI.ID files.
!! NOTE: vcoud and zcousq are in module m_zmelt.
          qxx=qibz_k
c          if(kx<=nqibz) qxx=qibz_k
c          if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
          ifvcoud = iopen('Vcoud.'//charnum5(kx),0,0,0)
          do
            read(ifvcoud) ngb0
            read(ifvcoud) qvv
            if(allocated(vcoud)) deallocate(vcoud)
            allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
            read(ifvcoud) vcoud
            read(ifvcoud) zcousq
            if(sum(abs(qvv-qxx))<1d-6) goto 1133
          enddo
          if(sum(abs(qvv-qxx))>1d-6) then
            write(6,*)'qvv =',qvv
            write(6,*)'qxx=',qxx,kx
            call rx( 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent')
          endif
 1133     continue
          if( ngb0/=ngb ) then  !sanity check
            write(6,*)' qxx ngb0 ngb=',qxx,ngb0,ngb
            call rx( 'hsfp0.m.f:ngb0/=ngb')
          endif
!! ppovlz is used in get_zmel
!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
!! ppovlz= O Zcousq
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
          call readppovl0(qibz_k,ngc,ppovl)
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
          ppovlz(nbloch+1:nbloch+ngc,:)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
c            write(6,*)'sumcheck ppovlz 00000 =',sum(abs(ppovlz(:,:)))
          deallocate(zcousq,ppovl)
!! === open WVR,WVI ===
          if(.not.exchange) then
            ifrcw  = iopen('WVR.'//charnum5(kx),0,-1,mrecl)
            ifrcwi = iopen('WVI.'//charnum5(kx),0,-1,mrecl)
          endif
          nrot=0
          do irot = 1,ngrp
c            if( kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if(kr==0) cycle   ! next irot
              qbz_kr= qbz (:,kr) 
c            else
c              kr=-99999         !for sanity check
c              qbz_kr= 0d0
c              if( wgt0(kx-nqibz,irot)==0d0 ) cycle ! next irot
c            endif
            nrot=nrot+1
          enddo  
!TIME1_01000 ":BeforDo1000"


!! === loop 1000 over rotations irot ===
          do 1000 irot = 1,ngrp
c            if( kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if(kr==0) cycle
              qbz_kr= qbz (:,kr) 
c            else
c              kr=-99999         !for sanity check
c              qbz_kr= 0d0
c              if( wgt0(kx-nqibz,irot)==0d0 ) cycle
c            endif

!TIME0_1010
!! no. occupied (core+valence) and unoccupied states at q-rk
            qk =  q - qbz_kr        
            call readeval(qk, isp, ekq) 
            ekc(nctot+1:nctot+nband) = ekq (1:nband)
            nt0 = nocc (ekc,ef,.true.,nstatetot)
            ddw= .5d0
c            if(GaussSmear()) ddw= 10d0
            ddw= 10d0
            efp= ef+ddw*esmr
            efm= ef-ddw*esmr
            nt0p = nocc (ekc,efp,.true.,nstatetot)
            nt0m = nocc (ekc,efm,.true.,nstatetot)
!! nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
!! nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W
            if(exchange) then
              nbmax = nt0p-nctot
            else
              nbmax = nband
              nbmxe = nocc (ekc,ebmx(2),.true.,nstatetot)-nctot
              nbmax  = min(nband,nbmx(2),nbmxe)
              if(initp) then
                write(6,*)' nbmax=',nbmax
                initp=.false.
              endif
            endif
c$$$!! ntqxx is number of bands for <i|sigma|j>.
c$$$            ntqxx = nocc (omega-eftrue,ebmx(1),.true.,ntq) 
c$$$!bug -ef is added jan2013
c$$$!previous version do not give wrong results, but inefficient.
c$$$            ntqxx = min(ntqxx, nbmx(1))
c$$$            if(ntqxx<nband) then
c$$$              do i=ntqxx,1,-1   !redudce ntqxx when band tops are degenerated. !sep2012
c$$$                if(omega(i+1)-omega(i)<1d-2) then
c$$$                  ntqxx=i-1
c$$$                else
c$$$                  exit
c$$$                endif
c$$$              enddo
c$$$            endif  
c$$$            nbandmx(ip)=ntqxx   !number of bands to be calculated Sep2012.

            ntqxx = nbandmx(ip) !mar2015
            if(debug) write(6,*)' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
            nstate = nctot + nbmax ! = nstate for the case of correlation

!! Get matrix element zmelt= rmelt + img*cmelt, defined in m_zmel.F---
c            if(debug) write(6,*)'zzBBB ppovlz =',sum(abs(ppovlz(:,:))),kx,irot
            if(allocated(zmel)) deallocate(zmel)
            if(allocated(zmeltt)) deallocate(zmeltt)
!TIME1_1010 "Beforeget_zmelt"
! this return zmeltt (for exchange), or zmel (for correlation)
!TIME0_1088
            call get_zmelt(exchange,q,kx,qibz_k,irot,qbz_kr,kr,isp,
     &       ngc,ngb,nbmax,ntqxx,nctot,ncc=0)
!TIME1_1088 "get_zmelt"

c$$$!! cccccccccc START: old version, instead of get_zmelt ccccccccccc
c$$$         call  readcphi(q, nlmto,isp, quu, cphikq)
c$$$         if(debug) write(6,*) ' sxcf: 2'
c$$$         do      it = 1,ntq
c$$$            itp        = itq(it)
c$$$            cphiq(1:nlmto,it) = cphikq(1:nlmto,itp)
c$$$            write(*,*)'svvvv ',it, itp, sum(cphiq(:,it))
c$$$         enddo
c$$$               write(*,*)'srrrrr 1c',sum(cphiq(:,1:ntq)),ntq
c$$$
c$$$         call dinv33(qbas,0,qbasinv,det)
c$$$         if(debug) write(6,*) ' sxcf: 1'
c$$$         if(allocated(expikt)) deallocate(expikt)
c$$$         allocate(expikt(natom))
c$$$ccccccccccccccccccccccccccccccccc
c$$$!!    rotate atomic positions invrot*R = R' + T
c$$$               invr  = invrot (irot,invg,ngrp)
c$$$               tr    = tiat(:,:,invr)
c$$$               iatomp= miat(:,invr)
c$$$               symope= symgg(:,:,irot)
c$$$               shtv  = matmul(symope,shtvg(:,invr))
c$$$!TIME1 "before ppbafp_v2"
c$$$!TIME0 
c$$$
c$$$!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
c$$$c               call ppbafp_v2 (irot,ngrp,isp,nsp,
c$$$c     i              il,in,im,nlnm, !w(i_mnl),
c$$$c     d              nl,nn_,nclass,nlnmx,
c$$$c     i              mdimx,lx,nx_,nxx_, !Bloch wave    
c$$$c     i              cgr, nl-1,  !rotated CG
c$$$c     i              ppbrd,      !radial integrals
c$$$c     o              ppb)
c$$$               ppb = ppbir(:,irot,isp)
c$$$!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
c$$$               qk =  q - qbz_kr !qbz(:,kr)
c$$$               call readcphi(qk, nlmto,isp, quu, cphikq)
c$$$!TIME1 "before expikt"
c$$$!TIME0
c$$$
c$$$!! ====================================================
c$$$!!     matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)>
c$$$!!     including the phase factor exp(ik.T)
c$$$!!     B(rot*k,r) = B(k,invrot*r)
c$$$!! ====================================================
c$$$!! phase factors expikt(ia) is for exp(ik.T(R))
c$$$               do ia = 1,natom
c$$$                  expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) 
c$$$               end do
c$$$!! matrix elements
c$$$!!   core
c$$$               nt   = nctot + nbmax ! = nstate for the case of correlation
c$$$               allocate( zzzmel(nbloch,nt,ntqxx))
c$$$               call psicb_v2  (icore,ncore,ntqxx,iclass,
c$$$     i              dreal(expikt(1:natom)),dimag(expikt(1:natom)),
c$$$     i              cphiq, 
c$$$     i              ppb,  
c$$$     i              nlnmv,nlnmc,mdim,
c$$$     i              imdim,iatomp,
c$$$     d              mdimx,nlmto,nbloch,nlnmx,nt,ntqxx,natom,nclass,
c$$$     d              nl,nnc,     
c$$$     o              zzzmel)
c$$$               if(debug)  write(6,*) ' sxcf_fal2sc: goto psi2bc1'
c$$$ccccccccc cccccccccccccccccccccccccccccccccccccccc
c$$$               write(*,*)'srrrrr 1',sum(cphikq(1:nlmto,1:ntq))
c$$$               write(*,*)'srrrrr 1',sum(cphiq(1:nlmto,1:ntq))
c$$$               write(*,*)'srrrrr 1',sum(ppb)
c$$$               write(*,*)'srrrrr 1',sum(expikt)
c$$$               write(*,*)'srrrrr 1',sum(zzzmel)
c$$$
c$$$!!   valence
c$$$               call psi2b_v2  (nbmax, ntqxx,iclass,
c$$$     i              dreal(expikt(1:natom)),dimag(expikt(1:natom)),
c$$$     i              cphikq,     !occ    q-rk nband
c$$$     i              cphiq,      !unocc  q    ntq
c$$$     i              ppb,  
c$$$     i              nlnmv,nlnmc,mdim,nctot,
c$$$     i              imdim,iatomp,
c$$$     d              mdimx,nlmto,nbloch,nlnmx, nband, nt,ntqxx,
c$$$     d              natom,nclass,
c$$$     o              zzzmel)
c$$$               if(verbose()>50) call timeshow("4 after psi2bc1")
c$$$c               if(debug2) then
c$$$                  write(6,"('sum of zmel abszmel=',4d23.16)") sum(zzzmel),sum(abs(zzzmel) )
c$$$c               end if
c$$$!TIME1 "bfore psi2b_v2"
c$$$!TIME0
c$$$!! -- IPW part.
c$$$               if(debug) write(6,*) ' sxcf_fal1: goto drvmelp2 xxx111'
c$$$               allocate(drealzzzmel(nbloch,nt,ntqxx),dimagzzzmel(nbloch,nt,ntqxx))
c$$$               drealzzzmel=dreal(zzzmel)
c$$$               dimagzzzmel=dimag(zzzmel)
c$$$               deallocate(zzzmel)
c$$$               allocate( rmelt(ngb, nctot+nbmax, ntqxx), ! nstate= nctot+nband
c$$$     &              cmelt(ngb, nctot+nbmax, ntqxx))
c$$$               call drvmelp2( q,             ntqxx, ! q in FBZ
c$$$     i              q-qbz_kr,  nbmax, ! q-rk
c$$$     i              qibz_k,  ! k in IBZ for mixed product basis. rk = symope(qibz_k) 
c$$$     i              isp,ginv,
c$$$     i              ngc,ngcmx, ngpmx,nband,itq, 
c$$$     i              symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
c$$$     i              drealzzzmel, dimagzzzmel, nbloch, nt,nctot,
c$$$     o              rmelt,cmelt)
c$$$               if(debug) write(6,*) ' sxcf_fal1: end of drvmelp2'
c$$$               deallocate(drealzzzmel,dimagzzzmel)
c$$$               if(verbose()>50) call timeshow("5 after drvmelp")
c$$$               if(nbcut/=0.and.(.not.exchange)) then
c$$$                  do it= nctot+1,nctot+min(nbcut,nbmax)
c$$$                     rmelt(:, it,:) =0d0
c$$$                     cmelt(:, it,:) =0d0
c$$$                  enddo
c$$$               endif
c$$$               write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)
c$$$
c$$$!TIME1 "after drvmelp2"
c$$$!! NOTE:=========================================
c$$$!! zmelt = rmelt(igb(qbz_kr), iocc(q), iunocc(q-qbz_kr)) + i* cmelt
c$$$!! iunocc: band index at target  q.
c$$$!! iocc:   band index at intermediate vector qk = q - qbz_kr 
c$$$!! igb: index of mixed product basis       at qbz_kr (or written as rk)
c$$$!!   igb=1,ngb 
c$$$!!   ngb=nbloch+ngc  ngb: # of mixed product basis
c$$$!!                   nbloch: # of product basis (within MTs)
c$$$!!                   ngc: # of IPW for the Screened Coulomb interaction.
c$$$!!                   igc is for given 
c$$$!! See readgeig in drvmelp2.
c$$$!! ==================================================
c$$$!! smbasis ---need to fix this
c$$$c$$$               if(smbasis()) then !
c$$$c$$$                  ntp0= ntqxx
c$$$c$$$                  nn= nnr(kx)
c$$$c$$$                  no= nor(kx)
c$$$c$$$                  allocate( pomat(nn,no) )
c$$$c$$$                  pomat= pomatr(1:nn,1:no,kx)
c$$$c$$$                  if( sum(abs(qibz_k-qrr(:,kx)))>1d-10 .and.kx <= nqibz ) then
c$$$c$$$                     call rx( 'qibz/= qrr')
c$$$c$$$                  endif
c$$$c$$$                  if(no /= ngb.and.kx <= nqibz) then
c$$$c$$$!!     A bit sloppy check only for kx<nqibz because qibze is not supplied...
c$$$c$$$                     write(6,"(' q  ngb  ',3d13.5,3i5)")  qibz_k,ngb
c$$$c$$$                     write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
c$$$c$$$                     call rx( 'x0kf_v2h: POmat err no/=ngb')
c$$$c$$$                  endif
c$$$c$$$                  if(timemix) call timeshow("xxx2222 k-cycle")
c$$$c$$$                  ngb = nn      ! Renew ngb !!!
c$$$c$$$                  allocate ( zmel  (nn, nctot+nbmax, ntp0) )
c$$$c$$$                  call matm( pomat, dcmplx(rmelt,cmelt), zmel, 
c$$$c$$$     &                 nn, no, (nctot+nbmax)*ntp0 )
c$$$c$$$                  deallocate(rmelt, cmelt)
c$$$c$$$                  allocate( rmelt(ngb, nctot+nbmax, ntp0), !ngb is reduced.
c$$$c$$$     &                 cmelt(ngb, nctot+nbmax, ntp0) )
c$$$c$$$                  rmelt = dreal(zmel)
c$$$c$$$                  cmelt = dimag(zmel)
c$$$c$$$                  deallocate(zmel,pomat)
c$$$c$$$               else
c$$$c$$$                  nn=ngb
c$$$c$$$                  no=ngb
c$$$c$$$               endif
c$$$               nn=ngb
c$$$               no=ngb
c$$$               if( oncew() ) then
c$$$                  write(6,"('ngb nn no=',3i6)") ngb,nn,no
c$$$               endif
c$$$               if(timemix) call timeshow("22222 k-cycle")
c$$$!! === End of zmelt ; we now have matrix element zmelt= rmelt + img* cmelt === 
c$$$               if(allocated(zzzmel))deallocate(zzzmel) !rmel,cmel)
c$$$               if(debug) write(6,*) ' sxcf: goto wtt'
c$$$               if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)
c$$$
c$$$!! === End of zmelt ; we now have matrix element zmelt= rmelt + img* cmelt === 
c$$$!! cccccccccc  END: old version, instead of get_zmelt ccccccccccc


!! --- wtt setcion ---
c$$$               if(bzcase()==2)then
c$$$                  if(kx<=nqibz) then
c$$$                     wtt = wk(kr)
c$$$                     if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
c$$$                  elseif(kx>nqibz) then !   wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
c$$$                     wtt= wgt0(kx-nqibz,irot)
c$$$                  endif
c$$$               else
c            if(kx<= nqibz) then !  wtx = 1d0
              wtt = wk(kr)
c            else                !  wtx = wgt0(kx-nqibz,irot) 
c              wtt = wk(1)*wgt0(kx-nqibz,irot)
c              if(abs(wk(1)-1d0/dble(nqbz))>1d-10) call rx( 'sxcf:wk(1) inconsistent')
c            endif
!!
            if(eibz4sig()) then
              wtt=wtt*nrkip(kx,irot,ip) 
            endif   

!!--------------------------------------------------------
!! --- exchange section ---
!!--------------------------------------------------------
            if(exchange) then   !At the bottom of this block, cycle do 1000 irot.
!! We use the matrix elements zmeltt. Now given by "call get_zmelt"
!!
c need to check following comments ----     
c     S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c     v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
c     
c     > z1p(j,n,t) = S[i=1,nbloch] <psi(q,t) | psi(q-rk,n) B(rk,i)> v(k)(i,j)
c     
c     --- screened exchange case
c     if(screen) then
c     allocate( zw (nblochpmx,nblochpmx))
c     ix = 1
c     ! write(*,*)(kx-2)*(nw_w+1)+ix
c     read(ifrcw,rec=((kx-2)*nw+ix)) zw  ! Readin W(0) - v         !sf 22May02
c     !nw is number of frequency points in general mesh: freq_r(nw), freq_r(1)=0
c     vcoul = vcoul + zw(1:ngb,1:ngb) !c  screen test
c     deallocate(zw)
c     endif
!TIME0_0130
              vc = vcoud(1)     ! save vcoud(1)
              if (kx == iqini) vcoud(1) = wklm(1)* fpi*sqrt(fpi) /wk(kx)
              allocate(z1r(ntqxx,ngb),z2r(ntqxx,ngb),w3pi(ntqxx,ntqxx))
              allocate(w3p(nctot+nbmax,ntqxx,ntqxx))
              do  it = 1, nctot+nbmax
                do  ivc = 1, ngb
                  do  itp = 1, ntqxx
                    z1r(itp,ivc) = zmeltt(it,itp,ivc) * vcoud(ivc)
                    z2r(itp,ivc) = zmeltt(it,itp,ivc)
                  enddo         ! ivc
                enddo           ! it
                call zgemm('N','C',ntqxx,ntqxx,ngb,(1d0,0d0),z1r,ntqxx,
     .           z2r,ntqxx,(0d0,0d0),w3pi,ntqxx)
C              call zprm('w3pi',w3p,ntqxx,ntqxx,ntqxx)
C             Faster, but harder to parallelize
!                        call zqsmpy(11,'N','C',ntqxx,ngb,z1r,ntqxx,z2r,ntqxx,
!     .                       (0d0,0d0),w3pi,ntqxx)
C             call zprm('w3pi',w3p,ntqxx,ntqxx,ntqxx)
                do  itp = 1, ntqxx
                  do itpp = 1, ntqxx
                    w3p(it,itp,itpp) = w3pi(itp,itpp)
                  enddo
                enddo
              enddo
              vcoud(1) = vc     !restore vcoud(1)
              deallocate(z1r,z2r,w3pi)
              if(verbose()>=30) call cputid2(' complete w3p',0)
              deallocate(zmeltt)
              if(debug) then
                do  it  = 1,nctot+nbmax; do  itp = 1,ntqxx
                  write(6,"(' w3p =',2i4,2d14.6)") it,itp,w3p(it,itp,itp)
                enddo;    enddo
              endif
!TIME1_0130 "end_of_w3p"

c$$$#else
c$$$!kino 2014-08-13  !$OMP parallel  private(vc) 
c$$$!kino 2014-08-13  !$OMP do 
c$$$                     do itp= 1,ntqxx
c$$$                        do it = 1,nctot+nbmax
c$$$                           do ivc=1,ngb
c$$$                              zmeltt(it,itp,ivc) =  sum( zmel(:,it,itp)* ppovlz(:,ivc) )
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$!kino 2014-08-13  !$OMP end do 
c$$$!kino 2014-08-13  !$OMP do 
c$$$                     do 992 itpp= 1,ntqxx
c$$$                        do 993 itp = 1,ntqxx
c$$$                           if(diagonly.and.(itpp/=itp)) cycle
c$$$!! sep2013t a test:c           if(itpp>ntqxxd .and.itp/=itpp) cycle
c$$$                           do 994 it  = 1,nctot+nbmax
c$$$                              w3p(it,itp,itpp) = 0d0
c$$$                              do ivc=1,ngb
c$$$                                 if(ivc==1.and.kx==iqini) then
c$$$                                    vc= wklm(1)* fpi*sqrt(fpi) /wk(kx)
c$$$c     write(6,*)'wklm(1) vc=',wklm(1),vc
c$$$                                 else
c$$$                                    vc= vcoud(ivc)
c$$$                                 endif
c$$$c     zmelt1 =  sum( zmel(:,it,itp)  *ppovlz(:,ivc) )
c$$$c     zmelt2 =  sum( zmel(:,it,itpp) *ppovlz(:,ivc) )
c$$$                                 w3p(it,itp,itpp) = w3p(it,itp,itpp)
c$$$     &                                + vc * zmeltt(it,itp,ivc)*dconjg(zmeltt(it,itpp,ivc))
c$$$                              enddo
c$$$ 994                       continue
c$$$ 993                    continue
c$$$ 992                 continue
c$$$!kino 2014-08-13  !$OMP end do 
c$$$!kino 2014-08-13  !$OMP end parallel 
c$$$#endif
!KINO                     write(*,*)'kino: w3p checksum=',sum(w3p)
c                     deallocate(zmeltt)
c$$$                  else
c$$$!kino 2014-08-13  !$OMP parallel do 
c$$$                     do itpp= 1,ntqxx
c$$$                        do itp = 1,ntqxx
c$$$                           if(diagonly.and.(itpp/=itp)) cycle
c$$$c sep2013t a test:c  if(itpp>ntqxxd .and.itp/=itpp) cycle
c$$$                           do it  = 1,nctot+nbmax
c$$$                              w3p(it,itp,itpp) =dcmplx(
c$$$     &                          sum ( dreal(z1p(:,it,itp))*rmelt(:,it,itpp)
c$$$     &                          +   dimag(z1p(:,it,itp))*cmelt(:,it,itpp) ) ,
c$$$     &                          sum ( dimag(z1p(:,it,itp))*rmelt(:,it,itpp) 
c$$$     &                          -   dreal(z1p(:,it,itp))*cmelt(:,it,itpp) ) )
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$!kino 2014-08-13  !$OMP end parallel do 
c$$$                     deallocate(z1p)
c$$$                  endif  
c                  deallocate(zmel)
c$$$!!-- Write the Spectrum function for exchange May. 2001
c$$$                  if(ifexsp/=0) then
c$$$                     do it  = 1, nctot+nbmax
c$$$                        do itp = 1,ntqxx
c$$$                           write(ifexsp,"(3i4, 3f12.4, ' ',d23.15,'  ',d23.15)")
c$$$     &                          ip,itp,it, qbz_kr, ekc(it), -wtt*dreal(w3p(it,itp,itp))
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$!TIME1 "end of write ifsexsp"

!TIME0_0180 
!! --- Correct weigts wfac for valence by esmr
              do it = nctot+1, nctot+nbmax
                wfac = wfacx(-1d99, ef, ekc(it), esmr) !gaussian
                w3p(it,1:ntqxx,1:ntqxx) = wfac * w3p(it,1:ntqxx,1:ntqxx)
              enddo

!! apr2015 correct weights for core-hole case
              if(corehole) then
              do it = 1, nctot
                w3p(it,1:ntqxx,1:ntqxx) = wcorehole(it,isp) * w3p(it,1:ntqxx,1:ntqxx)
              enddo
              endif

              do itpp=1,ntqxx
                do itp = 1,ntqxx !S[j=1,nbloch]  z1p(j,t,n) <B(rk,j) psi(q-rk,n) |psi(q,t')>
                  if(jobsw==5.and.(itpp/=itp)) cycle
c sep2013t a test:c  if(itpp>ntqxxd .and.itp/=itpp) cycle
                  zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
     &             - wtt * sum( w3p(:,itp,itpp) )
                enddo
              enddo
              deallocate( w3p)
c$$$                  if(.not.newaniso()) deallocate(vcoul)
!TIME1_0180 "enddo_zsec_wtt_sum"
              cycle             ! next irot do 1000 loop
            endif               ! end of if(exchange)
!! ============== End of exchange section ======================
            if(timemix) call timeshow("33333 k-cycle")
cc!TIME1 "end of exchange section"


!!----------------------------------------------------------
!!---  correlation section ---------------------------------
!!----------------------------------------------------------
!!  We use the matrix elements zmel, which is given by "call get_zmelt"
!!
!!================================================================
!! need to check the following notes.
!!     The correlated part of the self-energy:
!!     S[n=all] S[i,j=1,nbloch]
!!     <psi(q,t) |psi(q-rk,n) B(rk,i)>
!!     < [w'=0,inf] (1/pi) (w-e)/{(w-e)^2 + w'^2} Wc(k,iw')(i,j) >
!!     <B(rk,j) psi(q-rk,n) |psi(q,t)>
!!     e = e(q-rk,n), w' is real, Wc = W-v
!!================================================================
!! Get zwz0(omega=0, m, i, j), and zwz(i omega, m, i, j)
!! m intermediate state. zwz= \sum_I,J <i|m I> W_IJ(i omega) <J m|j>
!!
!! sum over both occupied and unoccupied states and multiply by weight
!     new from Jan2006! I think this should be OK.  ----------------------------
!     The output of sxcf_fal2 is  <i|Re[S](e_i)|j> ------------
!     Im-axis integral gives Hermitian part of S.
!     (Be careful as for the difference between
!     <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
!     ---because e_i is included.
!     The symmetrization (hermitian) procedure is inlucded in hqpe.sc.F
!     old befor Jan2006
!     &        wtt*.5d0*(   sum(zwzi(:,itp,itpp))+ !S_{ij}(e_i)
!     &        dconjg( sum(zwzi(:,itpp,itp)) )   ) !S_{ji}^*(e_j)= S_{ij}(e_j)
!-----------------------------------------------------------------------------
!! omega integlation along im axis.
!! zwzi(istate,itqxx1,itqxx2)=\int_ImAxis d\omega' zwz(omega',istate,itqxx1,itqxx2) 1/(omt-omega')
!! ,where omt=omegat is given in the following 1385-1386 loop.
!!


!! -------------------------------------------------------------------
!! Contribution to SEc(qt,w) from integration along the imaginary axis
!!     loop over w' = (1-x)/x, frequencies in Wc(k,w')
!!     {x} are gaussian-integration points between (0,1)
!!---------------------------------------------------------------------
!! Readin W(omega=0) and W(i*omega)
!! Then get zwz0 and zwz
!! zwz0 = (zmel*)*(W(*omega=0)   -v)*zmel
!! zwz =  (zmel*)*(W(i*omega(ix))-v)*zmel
!TIME0_0200 
            allocate( zwz0(        nstate,ntqxx,ntqxx))
            allocate( zwz (niw*npm,nstate,ntqxx,ntqxx)) 
            allocate( zw (nblochpmx,nblochpmx))
            ix = 1 + (0 - nw_i) !at omega=0 ! nw_i=0 (Time reversal) or nw_i =-nw
            read(ifrcw,rec=ix) zw ! direct access read Wc(0) = W(0) - v
            call matzwz2(2, zw(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb,   
     o       zwz0)                   
            do 1380 istate=1,nstate 
              zwz0(istate,1:ntqxx,1:ntqxx) = ! w(iw) + w(-iw) Hermitian part.
     &         (zwz0(istate,1:ntqxx,1:ntqxx)
     &         + dconjg(transpose(zwz0(istate,1:ntqxx,1:ntqxx))))/2d0 
 1380       continue
            do 1390 ix=1,niw    !niw is usually ~10 points.
              read(ifrcwi,rec=ix) zw ! direct access read Wc(i*omega)=W(i*omega)-v
              call matzwz2(2, zw(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb,   
     o         zwz(ix,1:nstate,1:ntqxx,1:ntqxx)) ! zwz = zmel*(W(0)-v)*zmel
              do 1395 istate=1,nstate
                zw(1:ntqxx,1:ntqxx)= zwz(ix,istate,1:ntqxx,1:ntqxx)
                zwz(ix,istate,1:ntqxx,1:ntqxx) = ! w(iw) + w(-iw)  Harmitian part
     &           ( zw(1:ntqxx,1:ntqxx)
     &           + dconjg(transpose(zw(1:ntqxx,1:ntqxx))) )/2d0 
                if(npm==2) then ! w(iw) - w(-iw) Anti Hermitian part
                  zwz(ix+niw,istate,1:ntqxx,1:ntqxx) = 
     &             ( zw(1:ntqxx,1:ntqxx)
     &             - dconjg(transpose(zw(1:ntqxx,1:ntqxx))) )/2d0/img
                endif
 1395         continue
 1390       continue
            deallocate(zw)
!TIME1_0200 "endofdo1390"
!! Integration along imag axis for zwz(omega) for given it,itp,itpp
!! itp  : left-hand end of expternal band index.
!! itpp : right-hand end of expternal band index.
!! it   : intermediate state of G.
!TIME0_0210
            allocate(zwzi(nstate,ntqxx,ntqxx)) 
            do 1400 itpp= 1,ntqxx
              do 1410 itp = 1,ntqxx
                if((jobsw==5).and.(itpp/=itp)) cycle
                if (jobsw==1.or.jobsw==4) then
                  omegat = ef
c               elseif (jobsw==2)               omegat=.5d0*(omega(itp)+omega(itpp))
                else
                  omegat = omega(itp)
                endif  
                do 1420  it = 1,nstate
                  we =.5d0*( omegat -ekc(it))
                  if(it <= nctot) then
                    esmrx = 0d0
                  else
                    esmrx = esmr
                  endif   
!! ua_auto may be recovered in future...
c     if(ua_auto) then
c     ratio = .5d0 *( abs(zwz(niw,it,itp,itp  )/zwz0(it,itp,itp  ))
c     &                     +abs(zwz(niw,it,itpp,itpp)/zwz0(it,itpp,itpp)) )
c     call gen_ua(ratio,niw,freqx,  expa_,ua_)
c     endif
!! Gaussian smearing. Integration along im axis. zwz(1:niw) and zwz0 are used.
                  zwzi(it,itp,itpp) = 
     &             wintzsg_npm (npm, zwz(1,it,itp,itpp), zwz0(it,itp,itpp) 
     &             ,freqx,wx,ua_,expa_,we,niw,esmrx)
c                    zwzi(it,itp,itpp) =  !rectangular smearing only for npm=1
c     &                          wintzav ( zwz(1,it,itp,itpp),zwz0(it,itp,itpp)
c     &                          ,freqx,wx,ua_,expa_,we,niw, esmrx)
 1420           continue
 1410         continue  
 1400       continue  
            deallocate(zwz0,zwz) !zwzs
            if(debug) print *,'zzzzzzzzz sum zwzi ',sum(abs(zwzi(:,:,:)))
!TIME1_0210 "endofdo1400"
!! Contribution to Sigma_{ij}(e_i)
            do  1500 itpp= 1,ntqxx
              do 1510 itp = 1,ntqxx
                if( jobsw==5.and.(itpp/=itp)) cycle
                zsec(itp,itpp,ip) = zsec(itp,itpp,ip) + wtt*sum(zwzi(:,itp,itpp)) 
 1510         continue
 1500       continue
            deallocate(zwzi)
            if(jobsw==4) goto 2002

!! -------------------------------------------------------------------
!!  Contribution to SEc(qt,w) from the poles of G (integral along real axis)
!!    Currently, jobsw =1,3,5 are allowed...
!!    The variable we means \omega_epsilon in Eq.(55) in PRB76,165106 (2007)
!! -------------------------------------------------------------------
!TIME0_0310 
            if(timemix) call timeshow("goto Sec pole part k-cycle")
            if(debug)  write(6,*)'GOTO contribution to SEc(qt,w) from the poles of G'
            if (.not.(jobsw == 1 .or. jobsw == 3.or.jobsw==5)) then
              call rx( 'sxcf_fal3_scz: jobsw /= 1 3 5')
            endif
!! Get index nwxi nwx nt_max. finish quickly. We can simplify this...
            call get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,
     i       nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,
     o       nwxi,nwx,nt_max)
!! assemble small arrays first.
            allocate(we_(nt_max,ntqxx),wfac_(nt_max,ntqxx),ixss(nt_max,ntqxx),ititpskip(nt_max,ntqxx),iirx(ntqxx))
            call weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,
     i       ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt,
     o       we_,wfac_,ixss,ititpskip,iirx)

!! We need zw3, the Hermitian part, because we need only hermitean part of Sigma_nn'
!! This can be large array; nwx-nwxi+1 \sim 400 or so...
            allocate( zw3(ngb,ngb,nwxi:nwx)) 
            allocate( zw(nblochpmx,nblochpmx))
            do ix = nwxi,nwx
              nrec = ix-nw_i+1  !freq_r(ix is in nw_i:nx)
              read(ifrcw,rec=nrec) zw ! direct access Wc(omega) = W(omega) - v
              if(hermitianW) then
                zw3(:,:,ix)=(zw(1:ngb,1:ngb)+transpose(dconjg(zw(1:ngb,1:ngb))))/2d0
              else
                zw3(:,:,ix)=zw(1:ngb,1:ngb)
              endif
            enddo
            deallocate(zw)
!! rearrange index of zmel
            allocate(zmel1(ngb)) 
            if(jobsw==3) then
              allocate(zmel1_(ntqxx,ngb,nstate))
              do itpp= 1,ntqxx
                do it  = 1,nstate
                  zmel1_(itpp,1:ngb,it) = zmel(1:ngb,it,itpp)
                enddo
              enddo
            endif
!! jobsw==3
            if( jobsw==3) then
              allocate(zwz44(3,ntqxx),zwz4(ntqxx,3))
              do itp=1,ntqxx
                do it=1,nt_max
                  if(ititpskip(it,itp)) cycle
                  we =  we_(it,itp) 
                  ixs=  ixss(it,itp)
                  zmel1(:)=dconjg(zmel(:,it,itp))
                  zwz4=0d0
                  do ix0=1,3
                    ix=ixs+ix0-2
                    do igb2=1,ngb 
! !                     **** most time consuming part ******
                      zz2=sum(zmel1(1:ngb)*zw3(1:ngb,igb2, iirx(itp)*ix)  )
                      call zaxpy(ntqxx,zz2,zmel1_(1,igb2,it),1,zwz4(1,ix0),1)
                    enddo
                  enddo  
                  zwz44 = transpose(zwz4)
                  do itpp=1,ntqxx
                    if(npm==1) then
                      zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
     .                 + wfac_(it,itp) * alagr3z2(we,freq_r(ixs-1),zwz44(1,itpp),itp==itpp ) !mar015 ,itp,itpp)
                    else
                      zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
     .                 + wfac_(it,itp) * alagr3z(we,freq_r(ixs-1),zwz44(1,itpp))
                    endif
                  enddo
                enddo
              enddo
              deallocate(zwz44,zwz4)
            endif

!! jobsw=1,5 Sigma are calculated.
            if( jobsw==1.or.jobsw==5) then
              do itp=1,ntqxx
                do it=1,nt_max
                  if(ititpskip(it,itp)) cycle
                  we =  we_(it,itp) 
                  ixs=  ixss(it,itp)
                  zmel1(:)=dconjg(zmel(:,it,itp))
                  zwz3=0d0
                  do ix0=1,3
                    ix=ixs+ix0-2
!!               **** most time consuming part for jobsw=1 ******
!!               To reduce computational time, confusing treatment only uses lower half of zw3 (zw3 is Hermitan)
!!               Clean up needed.

!! zwz3 contains <itp| it I> wz3_IJ(we)  <J it| itp>
!!    when zw3 is hermitian.
                    if(hermitianW) then
                      do igb2=2,ngb 
                      zz2 = sum(zmel1(1:igb2-1)*zw3(1:igb2-1,igb2,iirx(itp)*ix)  ) +
     &                 .5d0* zmel1(igb2)*zw3(igb2,igb2,iirx(itp)*ix)
                      zwz3(ix0) = zwz3(ix0)+zz2*zmel(igb2,it,itp)
                      enddo       !igb2
                      zwz3(ix0) = 2d0*dreal(zwz3(ix0))+ !I think 2d0 is from upper half.
     &                zmel1(1)*zw3(1,1, iirx(itp)*ix)*zmel(1,it,itp)
!!    when zw3 is not need to be hermitian case. This gives life time
                    else  
                      zwz3(ix0) = sum( matmul(zmel1(1:ngb), zw3(1:ngb,1:ngb,iirx(itp)*ix))*zmel(1:ngb,it,itp) )
                    endif  
                  enddo
                  if(npm==1) then
                    zsec(itp,itp,ip) = zsec(itp,itp,ip) 
     .               + wfac_(it,itp)*alagr3z2(we,freq_r(ixs-1),zwz3,.true.)
                  else
                    zsec(itp,itp,ip) = zsec(itp,itp,ip) 
     .               + wfac_(it,itp)*alagr3z(we,freq_r(ixs-1),zwz3)
                  endif
                enddo
              enddo
            endif               
!TIME1_0310 "EndReCorrelation"
c           goto 2012

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c$$$
c$$$
c$$$
c$$$ccccccccc old code cccccccccccccccccccccccc
c$$$            if(timemix) call timeshow("55555 k-cycle")
c$$$            if(debug) write(*,'(a,5i6)')'kino: ntqxx,itini,itend,ngb=',ntqxx,itini,itend,ngb
c$$$c$$$               if(test_symmetric_W().and.npm==2) then
c$$$c$$$                  if(onceww(4)) write(6,*)' test_symmetric_W()=',test_symmetric_W(),nwxi,nwx
c$$$c$$$                  allocate(zw3x(ngb,ngb))
c$$$c$$$                  do ix= 1,min(abs(nwxi),nwx)
c$$$c$$$                     zw3x = 0.5d0* (zw3(:,:,ix) + zw3(:,:,-ix))
c$$$c$$$                     zw3(:,:, ix)=zw3x
c$$$c$$$                     zw3(:,:,-ix)=zw3x
c$$$c$$$                  enddo
c$$$c$$$                  deallocate(zw3x)
c$$$c$$$               endif
c$$$!TIME1 "before 2001"
c$$$!TIME0 
c$$$            allocate(zwz44(3,ntqxx),zwz4(ntqxx,3))
c$$$            do 2001 itp = 1,ntqxx ! loop over states (q-k,n)
c$$$              omg = omega(itp)
c$$$              if (omg >= ef) then
c$$$                itini= nt0m+1
c$$$                itend= nt_max
c$$$                iii=  1
c$$$              else
c$$$                itini= 1
c$$$                itend= nt0p
c$$$                iii= -1
c$$$              endif
c$$$              do 2011 it = itini,itend ! nt0p corresponds to efp
c$$$                esmrx = esmr
c$$$                if(it<=nctot) esmrx = 0d0
c$$$                wfac = wfacx2(omg,ef, ekc(it),esmrx)
c$$$                if(wfac<wfaccut) cycle ! next it
c$$$                we = .5d0* abs( omg-weavx2(omg,ef, ekc(it),esmr) ) !Gaussian smearing
c$$$                if(it<=nctot .and.wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
c$$$c$$$                    Rectangular smearing
c$$$c$$$                        if( wfac==0d0) cycle ! next it
c$$$c$$$                        if( omg >= ef) we = 0.5d0* abs( max(omg-ekc(it), 0d0) )
c$$$c$$$                        if( omg <  ef) we = 0.5d0* abs( min(omg-ekc(it), 0d0) )
c$$$c$$$                        if( it<=nctot) then !faleev
c$$$c$$$                           if(wfac/=0) call rx( "sxcf:  it<=nctot.and.wfac/=0")
c$$$c$$$                        endif
c$$$c$$$                     endif
c$$$                if(debug) write(6,"( ' xxx1',10d13.6)") omg,ef, ekc(it),wfac
c$$$                wfac= iii* wfac*wtt
c$$$                do iwp = 1,nw 
c$$$                  ixs=iwp
c$$$                  if(freq_r(iwp)>we) exit
c$$$                enddo  
c$$$                if(nw_i==0) then
c$$$                  if(ixs+1>nwx) call rx( ' sxcf: ixs+1>nwx xxx2')
c$$$                else
c$$$                  if(omg >=ef .and. ixs+1> nwx ) then
c$$$                    write(6,*)'ixs+1 nwx=',ixs+1,nwx
c$$$                    call rx( ' sxcf: ixs+1>nwx yyy2a')
c$$$                  endif
c$$$                  if(omg < ef .and. abs(ixs+1)> abs(nwxi) ) then
c$$$                    write(6,*)'ixs+1 nwxi=',ixs+1,nwxi
c$$$                    call rx( ' sxcf: ixs-1<nwi yyy2b')
c$$$                  endif
c$$$                endif
c$$$                iir = 1
c$$$                if(omg < ef .and. nw_i/=0) iir = -1
c$$$                zmel1(:)=dconjg(zmel(:,it,itp))
c$$$
c$$$                if (jobsw == 1.or.jobsw==5) then
c$$$                  zwz3=(0d0,0d0)
c$$$!kino 2014-08-13  !$OMP parallel do private(ix,zz2)
c$$$                  do 2014 ix0=1,3
c$$$                    ix=ixs+ix0-2
c$$$                    do igb2=2,ngb 
c$$$! !**** most time consuming part for jobsw=1 ******
c$$$                      zz2=sum(zmel1(1:igb2-1)*zw3(1:igb2-1,igb2,iir*ix)  ) +
c$$$     &                 .5d0* zmel1(igb2)*zw3(igb2,igb2,iir*ix)
c$$$                      zwz3(ix0)=zwz3(ix0)+zz2*zmel(igb2,it,itp)
c$$$                    enddo       !igb2
c$$$                    zwz3(ix0)=2d0*dreal(zwz3(ix0))+
c$$$     &               zmel1(1)*zw3(1,1, iir*ix)*zmel(1,it,itp)
c$$$ 2014             continue      !ix
c$$$!kino 2014-08-13  !$OMP end parallel do 
c$$$                  if(npm==1) then
c$$$                    zsec(itp,itp,ip) = zsec(itp,itp,ip) 
c$$$     .               + wfac*alagr3z2(we,freq_r(ixs-1),zwz3,itp,itp)
c$$$                  else
c$$$                    zsec(itp,itp,ip) = zsec(itp,itp,ip) 
c$$$     .               + wfac*alagr3z(we,freq_r(ixs-1),zwz3)
c$$$                  endif
c$$$!!  this contribution to zsec_nn is real (hermitean)
c$$$
c$$$                elseif(jobsw == 3) then
c$$$                  zwz4=(0d0,0d0)
c$$$!$OMP parallel private(ix,zz2)
c$$$                  do 2015 ix0=1,3
c$$$                    ix=ixs+ix0-2
c$$$!$OMP do reduction(+:zwz4)
c$$$!! Next zaxpy is most time consuming part for jobsw=3.****
c$$$!! I think we can speed up this section...
c$$$                    do igb2=1,ngb 
c$$$                      zz2=sum(zmel1(1:ngb)*zw3(1:ngb,igb2, iir*ix)  )
c$$$                      call zaxpy(ntqxx,zz2,zmel1_(1,igb2,it),1,zwz4(1,ix0),1)
c$$$                    enddo
c$$$ 2015             continue      !ix0
c$$$!$OMP end parallel
c$$$                  zwz44 = transpose(zwz4)
c$$$                  do itpp=1,ntqxx
c$$$                    if(jobsw==5.and.(itpp/=itp)) cycle
c$$$                    if(npm==1) then
c$$$                      zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
c$$$     .                 + wfac*alagr3z2(we,freq_r(ixs-1),zwz44(1,itpp),itp,itpp)
c$$$                    else
c$$$                      zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
c$$$     .                 + wfac*alagr3z(we,freq_r(ixs-1),zwz44(1,itpp))
c$$$                    endif
c$$$                  enddo         !itpp
c$$$                endif           ! inner jobsw=1 or 3
c$$$!!     this contribution to zsec_nn' is not hermitean because W(e_n)
c$$$!!     and must be made hermitean when zsec will be written on disc
c$$$ 2011         continue
c$$$ 2001       continue            !itp
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc

 2012       continue
            deallocate(we_,wfac_,ixss,ititpskip,iirx)
 2002       continue
            deallocate(zw3, zmel, zmel1)
            if(allocated(zmel1_)) deallocate(zmel1_)
            if(verbose()>50) call timeshow("11after alagr3z iw,itp,it cycles")
            if(debug) then
              write(6,*)' end of do 2001 '
              do itp = 1,ntq
                write(6,'(" zsec=",i3,2d15.7)') itp,zsec(itp,itp,ip)
              enddo
            endif
 1000     continue              ! end do irot
          ifvcoud =iclose('Vcoud.'//charnum5(kx))
          if(.not.exchange) then
            ifrcw  = iclose('WVR.'//charnum5(kx))
            ifrcwi = iclose('WVI.'//charnum5(kx))
          endif   
          deallocate(ppovlz)
 1100   continue                ! end of kx-loop
        ifvcoud =iclose('Vcoud')
        if(irot==1) write(6,"('  sum(abs(zsec))=',d23.15)") sum(abs(zsec))
        if (allocated(vcoul))deallocate(vcoul)
 1001 continue                  ! end do ip
c      if(allocated(freq_r))deallocate(freq_r)
c      if (allocated(expikt))deallocate(expikt)
c!TIME1_0000 "end of sxcf_fal3_scz"
      end subroutine sxcf_fal3_scz


      subroutine weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,
     i ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt,
     o we_,wfac_,ixss,ititpskip,iirx)
!! generate required data set for main part of real part integration.
      implicit none
      integer,intent(in):: ntqxx,nctot,nw_i,nw,nt0m,nwx,nwxi,nt_max
      real(8),intent(in)::ef,omega(ntqxx),ekc(ntqxx),freq_r(nw_i:nw),esmr,wfaccut,wtt
      real(8),intent(out):: we_(nt_max,ntqxx),wfac_(nt_max,ntqxx)
      integer,intent(out) :: ixss(nt_max,ntqxx),iirx(ntqxx)
      logical,intent(out) :: ititpskip(nt_max,ntqxx)
      integer:: itini,iii,it,itend,wp,ixs,itp,iwp,nt0p
      real(8):: omg,esmrx,wfacx2,we,wfac,weavx2
      ititpskip=.false.
      do itp = 1,ntqxx          !this loop should finish in a second
        omg = omega(itp)
!  jobsw==2  
!            if (jobsw==2)  omg=.5d0*(omega(itp)+omega(itpp))
        iirx(itp) = 1
        if( omg < ef .and. nw_i/=0) iirx(itp) = -1
        if (omg >= ef) then
          itini= nt0m+1
          itend= nt_max
          iii=  1
        else
          itini= 1
          itend= nt0p
          iii= -1
        endif
        ititpskip(:itini-1,itp)=.true.
        ititpskip(itend+1:,itp)=.true.
        do it = itini,itend     ! nt0p corresponds to efp
          esmrx = esmr
          if(it<=nctot) esmrx = 0d0
          wfac_(it,itp) = wfacx2(omg,ef, ekc(it),esmrx)
          wfac = wfac_(it,itp)
          if(wfac<wfaccut) then
            ititpskip(it,itp)=.true.
            cycle 
          endif  
          wfac_(it,itp)=  wfac_(it,itp)*wtt*iii
!   Gaussian smearing we_= \bar{\omega_\epsilon} in sentences next to Eq.58 in PRB76,165106 (2007)
!   wfac_ = $w$ weight (smeared thus truncated by ef). See the sentences.
          we_(it,itp) = .5d0* abs( omg-weavx2(omg,ef, ekc(it),esmr) ) 
          we= we_(it,itp) 
          if(it<=nctot .and.wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
          do iwp = 1,nw 
            ixs = iwp
            if(freq_r(iwp)>we) exit
          enddo
          ixss(it,itp) = ixs
          if(nw_i==0) then
            if(ixs+1>nwx) call rx( ' sxcf: ixs+1>nwx xxx2')
          else
            if(omg >=ef .and. ixs+1> nwx ) then
              write(6,*)'ixs+1 nwx=',ixs+1,nwx
              call rx( ' sxcf: ixs+1>nwx yyy2a')
            endif
            if(omg < ef .and. abs(ixs+1)> abs(nwxi) ) then
              write(6,*)'ixs+1 nwxi=',ixs+1,nwxi
              call rx( ' sxcf: ixs-1<nwi yyy2b')
            endif
          endif
        enddo     
      enddo
      end subroutine weightset4intreal
      end module m_sxcfsc
!! --------------------------------------------------------------
      subroutine get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,
     i nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,
     o nwxi,nwx,nt_max)
!> Determine indexes of a range for calculation.
!! It is better to clean this up...
      implicit none
      integer,intent(in) :: nctot,nw_i,nw,nstate,nt0p,nt0m,ntq,
     & nband,ntqxx
      real(8),intent(in):: omega(ntq),esmr,ef,ekc(nctot+nband),wfaccut,
     & freq_r(nw_i:nw)
      integer,intent(out) :: nt_max,nwxi,nwx

      integer:: itp,it,itini,itend,iwp,ixs,ixsmin,ixsmx,verbose
      real(8):: omg,wfac,wfacx2,we,weavx2,esmrx,wexx
      logical::debug
!!     maximum ixs reqired.
      ixsmx =0
      ixsmin=0
      do 301 itp = 1,ntqxx
        omg  = omega(itp) 
        if (omg < ef) then
          itini= 1
          itend= nt0p
        else
          itini= nt0m+1
          itend= nstate
        endif
        do 311 it=itini,itend
          esmrx = esmr
          if(it<=nctot) esmrx = 0d0
          wfac = wfacx2(omg,ef, ekc(it),esmrx)
          if(wfac<wfaccut) cycle !Gaussian case
          we = .5d0*(weavx2(omg,ef,ekc(it),esmr)-omg)
cc Gaussian=F case  keep here just as a memo
c           if(wfac==0d0) cycle ! next it
c           if(omg>=ef) we = max( .5d0*(omg-ekc(it)), 0d0) ! positive
c           if(omg< ef) we = min( .5d0*(omg-ekc(it)), 0d0) ! negative
          if(it<=nctot) then
            if(wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
          endif
          do iwp = 1,nw
            ixs=iwp
            if(freq_r(iwp)>abs(we)) exit
          enddo
c     This change is because G(omega-omg') W(omg') !may2006
c     if(ixs>ixsmx  .and. omg<=ef ) ixsmx  = ixs
c     if(ixs>ixsmin .and. omg> ef ) ixsmin = ixs
          if(ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
          if(ixs>ixsmin .and. omg< ef ) ixsmin = ixs
          wexx  = we
          if(ixs+1 > nw) then
            write (*,*) ' nw_i ixsmin',nw_i, ixsmin
            write (*,*) ' wexx ',wexx
            write (*,*) ' omg ekc(it) ef ', omg,ekc(it),ef
            call rx( ' sxcf 222: |w-e| out of range')
          endif
 311    continue
 301  continue                  !end of SEc w and qt -loop
!!
      if(nw_i==0) then          !time reversal
        nwxi = 0
        nwx  = max(ixsmx+1,ixsmin+1)
      else                      !no time revarsal
        nwxi = -ixsmin-1
        nwx  =  ixsmx+1
      endif
      if (nwx > nw   ) then
        call rx( ' sxcf_fal3_sc nwx check : |w-e| > max(w)')
      endif
      if (nwxi < nw_i) then
        call rx( ' sxcf_fal3_sc nwxi check: |w-e| > max(w)')
      endif
      if(debug) write(6,*)'nw, nwx=',nw,nwx
      if(verbose()>50)call timeshow("10before alagr3z iw,itp,it ")
!!  Find nt_max 
      nt_max=nt0p               !initial nt_max
      do 401 itp = 1,ntqxx
        omg     = omega(itp)
        if (omg > ef) then
          do  it = nt0m+1,nstate ! nt0m corresponds to efm
            wfac = wfacx2 (ef,omg, ekc(it),esmr)
c                        if( (GaussSmear().and.wfac>wfaccut)
c     &                       .or.(.not.GaussSmear().and.wfac/=0d0)) then
            if(wfac>wfaccut) then
              if (it > nt_max) nt_max=it ! nt_max is  unocc. state
            endif               ! that ekc(it>nt_max)-omega > 0
          enddo                 ! so it > nt_max does not contribute to omega pole integral
        endif
 401  continue                  !end of  w and qt -loop
      end subroutine get_nwx
