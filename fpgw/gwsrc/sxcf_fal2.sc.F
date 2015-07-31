      module m_sxcf
      contains
      subroutine sxcf_fal3_scz(kount,qip,itq,ntq,ef,esmr, !tiat,miat,
     i nsp,isp,             
     i qbas,ginv,
     i qibz,qbz,wk,nstbz,irkip,nrkip, 
     i freq_r,nw_i,nw, freqx,wx,dw,         
     i ecore,
     i nlmto,nqibz,nqbz,nctot,
     i nl,nnc,nclass,natom,
     i nbloch,ngrp,niw,nq,      !nlnmx,mdimx,
     i nblochpmx ,ngpmx,ngcmx, 
     i wgt0,nq0i,q0i,symgg, alat, nband, ifvcfpout,
     i exchange,screen,cohtest, ifexsp,
     i nbmx,ebmx,           
     i wklm,lxklm,          
     i invg,il,in,im,nn_,lx,nx_,nxx_, !ppbrd, !,nlnm
     i eftrue,
!!
     i jobsw,                   != iSigma_en 
     o zsec,coh,nbandmx)
c     o exx,zsecd)

      use m_mpi
      use m_readqg,only:readqg0
      use m_readeigen,only: readeval,readcphi
      use keyvalue,only: getkeyvalue
      use m_zmel,only: get_zmelt,ppovlz,  zmel,zmeltt 
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
!! \endverbatim
!!
!!
!! Calculate the exchange part and the correlated part of self-energy.
!! T.Kotani started development after the analysis of F.Aryasetiawan's LMTO-ASA-GW.
!! We still use some of his ideas in this code.
!!
!! See paper   
!! [1]T. Kotani and M. van Schilfgaarde, ??Quasiparticle self-consistent GW method: 
!!     A basis for the independent-particle approximation,?? Phys. Rev. B, vol. 76, no. 16, p. 165106[24pages], Oct. 2007.
!! [2]T. Kotani, ??Quasiparticle Self-Consistent GW Method Based on the Augmented Plane-Wave 
!!    and Muffin-Tin Orbital Method,?? J. Phys. Soc. Jpn., vol. 83, no. 9, p. 094711 [11 Pages], Sep. 2014.
!!
!! =======================================================================
!! Omega integral for SEc
!! The integral path is deformed along the imaginary-axis, but together with contribution of poles.
!! See Fig.1 and around in Ref.[1].
!!
!! ---Integration along imaginary axis.---
!!( Current version for it, wintzsg_npm, do not assume time-reversal when npm=2.)
!! Integration along the imaginary axis: -----------------
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
      integer:: dummy4doxygen

! input variables
      logical, intent(in) :: exchange,screen,cohtest
      integer, intent(in) :: ntq,natom,nqbz,nqibz,ngrp,nq,niw
      integer, intent(in) :: nband,nlmto,nq0i,nctot,isp,nsp !,mdim(*) !,nlnmx
      integer, intent(in) :: ifvcfpout,nbloch,nclass,nl,nnc, nblochpmx
      integer, intent(in) :: itq(ntq) !,nstar(nqibz) !miat(natom,ngrp),mdimx,
      integer, intent(in) :: irkip(nqibz,ngrp,nq),nrkip(nqibz,ngrp,nq)
      integer, intent(in) :: kount(nqibz,nq),ngpmx,ngcmx,ifexsp,jobsw
      integer, intent(in) :: nbmx(2) !,nlnmv(*),nlnmc(*)!,iclass(*),icore(*)
      integer, intent(in) :: nstbz(nqbz) !,nomx !,nkpo,nnmx,imdim(*)ncore(*),
      integer, intent(in) :: lxklm,invg(ngrp) !nnr(nkpo),nor(nkpo),
      integer, intent(in) :: il(*),in(*),im(*),nn_,lx(*),nx_(*),nxx_ !,nlnm(*)
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
      integer :: iatomp(natom), nt0p, nt0m,nstate , nbmax, ntqxx
      integer :: nt,nw,ixs,iw,ivc,ifvcoud,ngb0
      integer :: nprecx,mrecl,ifwd,nrot,nwp,nw_i,ierr 
      integer :: nstatetot,iqini,iqend, ngb,ngc !nbcut,
      integer :: invr,nbmxe,ia,nn,ntp0,no,itpp,nrec,npm,itini,itend
      integer :: iwp,nwxi,nwx,iir, igb1,igb2,ix0,iii

      real(8) :: tpi, ekc(nctot+nband),ekq(nband), det, q(3),ua_
      real(8) :: expa_(niw), qxx(3), tr(3,natom), symope(3,3),shtv(3)
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
      complex(8),allocatable :: zw_zmel(:,:), zmel_zw_zmel(:), work(:)
      complex(8),allocatable :: zwz(:,:,:,:), zwz0(:,:,:),zwzi(:,:,:)
      complex(8),allocatable :: zwix(:,:),zwzix(:,:,:),zmel1(:) !,expikt(:)
      complex(8), allocatable :: zmel1_(:,:,:), zw3(:,:,:),zw3x(:,:)
      complex(8), allocatable :: zwz4(:,:),zwz44(:,:),pomat(:,:), zwzs(:)
      complex(8),allocatable :: ppovl(:,:),zcousq(:,:)
      complex(8),allocatable :: z1r(:,:),z2r(:,:),w3pi(:,:)
      complex(8) :: zmel_vc_zmel

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
c      real(8),allocatable:: ppb(:)
c      allocate( ppb(nlnmx*nlnmx*mdimx*nclass))

c      real(8)::exxq
      complex(8) :: ZDOTC, ZDOTU


c--------------------------------------------------------------------
!TIME0 
!TIME0 
!TIME0 
c      write(6,*)'sxcf_fal3_scz'
      timemix=.false.
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      debug=.false.
      if(verbose()>=90) debug=.true.

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
!TIME1 "before 1001"
      do 1001 ip = 1,nq         
        if(sum(irkip(:,:,ip))==0) cycle ! next ip

        write (6,*) ip,'  out of ',nq,'  k-points(extrnal q) '
!TIME0 
        q(1:3)= qip(1:3,ip)
        call readeval(q,isp,eq)
        do i  = 1,ntq
          omega(i) = eq(itq(i)) 
        enddo

!! we only consider bzcase()==1
        if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
        if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')

!TIME1 "before 1100"
!! NOTE total number of 
!!    kx loop(do 1100) and irot loop (do 1000) makes all the k mesh points.
!!    When iqini=1 (Gamma point), we use effective W(q=0) defined in the paper.
        iqini=1
        iqend=nqibz             !no sum for offset-Gamma points.

        do 1100 kx = iqini,iqend 
          if(sum(irkip(kx,:,ip))==0) cycle ! next kx
!TIME0 
          write(6,*) ' ### do 1100 start kx=',kx,' from ',iqini,' through', iqend
          if( kx <= nqibz ) then
            qibz_k= qibz(:,kx)
          else
            qibz_k= 0d0
          endif
          if(verbose()>=40)  write(6,*) ' sxcf_fal2_scz: loop 1100 kx=',kx
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
          if(kx<=nqibz) qxx=qibz_k
          if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
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
!! used in get_zmel
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
            if( kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if(kr==0) cycle   ! next irot
              qbz_kr= qbz (:,kr) 
            else
              kr=-99999         !for sanity check
              qbz_kr= 0d0
              if( wgt0(kx-nqibz,irot)==0d0 ) cycle ! next irot
            endif
            nrot=nrot+1
          enddo  
!TIME1 "before 1000"


!! === loop 1000 over rotations irot ===
          do 1000 irot = 1,ngrp
!TIME0

            if( kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if(kr==0) cycle
              qbz_kr= qbz (:,kr) 
            else
              kr=-99999         !for sanity check
              qbz_kr= 0d0
              if( wgt0(kx-nqibz,irot)==0d0 ) cycle
            endif

!! no. occupied (core+valence) and unoccupied states at q-rk
            qk =  q - qbz_kr        
            call readeval(qk, isp, ekq) 
            ekc(nctot+1:nctot+nband) = ekq (1:nband)
            nt0 = nocc (ekc,ef,.true.,nstatetot)
            ddw= .5d0
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

            ntqxx = nbandmx(ip) !mar2015
            if(debug) write(6,*)' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
            nstate = nctot + nbmax ! = nstate for the case of correlation

            call MPI__sxcf_rankdividerPT( ntqxx, nctot, nbmax )

!! Get matrix element zmelt= rmelt + img*cmelt, defined in m_zmel.F---
            if(allocated(zmel)) deallocate(zmel)
            if(allocated(zmeltt)) deallocate(zmeltt)

      timer_s(6) = MPI_Wtime()
! this return zmeltt (for exchange), or zmel (for correlation)
            call get_zmelt(exchange,q,kx,qibz_k,irot,qbz_kr,kr,isp,
     &           ngc,ngb,nbmax,ntqxx)
      timer_e(6) = MPI_Wtime() + timer_e(6) - timer_s(6)


!! --- wtt setcion ---
            if(kx<= nqibz) then !  wtx = 1d0
              wtt = wk(kr)
            else                 !  wtx = wgt0(kx-nqibz,irot) 
              wtt = wk(1)*wgt0(kx-nqibz,irot)
              if(abs(wk(1)-1d0/dble(nqbz))>1d-10) call rx( 'sxcf:wk(1) inconsistent')
            endif

            if(eibz4sig()) then
              wtt=wtt*nrkip(kx,irot,ip) 
            endif   

!!--------------------------------------------------------
!! --- exchange section ---
!!--------------------------------------------------------
!TIME0
               if(exchange) then !At the bottom of this block, cycle do 1000 irot.
!! We use the matrix elements zmeltt. Now given by "call get_zmelt"
!!
      timer_s(7) = MPI_Wtime()
                  allocate( work(ntqxx) )


                  vc = vcoud(1) ! save vcoud(1)
                  if (kx == iqini) vcoud(1) = wklm(1)* fpi*sqrt(fpi) /wk(kx)

                  do itpp=mpi__Ps, mpi__Pe ! P parallelization
c$$$              do itpp=1,ntqxx

                     work(1:ntqxx) = 0.0d0
                     do itp = 1,ntqxx
                        if(jobsw==5.and.(itpp/=itp)) cycle

                        do it=mpi__Ts, mpi__Te ! T parallelization
c$$$                    do it=1, nctot+nbmax
                           if( it>=nctot+1 ) then
                              wfac = wtt * wfacx(-1d99, ef, ekc(it), esmr) !gaussian
                           else
                              wfac = wtt
                           end if

#ifdef USE_GEMM_FOR_SUM
                           zmel_vc_zmel = ZDOTC( ngb, zmeltt(1:ngb,it,itpp), 1, vcoud(1:ngb) * zmeltt(1:ngb,it,itp), 1 )
#else
! not USE_GEMM_FOR_SUM
                           zmel_vc_zmel = sum( dconjg(zmeltt(1:ngb,it,itpp))
     &                          *              vcoud(1:ngb)
     &                          *             zmeltt(1:ngb,it,itp) )
#endif

                           work(itp) = work(itp) + wfac * zmel_vc_zmel
                        end do ! it
                     end do ! itp

                     call MPI__AllreduceSumT(work,ntqxx)

                     do itp = 1,ntqxx
                        zsec(itp,itpp,ip) = zsec(itp,itpp,ip) - work(itp)
                     end do ! itp
                  end do ! itpp

                  deallocate(work)
                  deallocate(zmeltt)

                  vcoud(1) = vc !restore vcoud(1)
      timer_e(7) = MPI_Wtime() + timer_e(7) - timer_s(7)
                  cycle         ! next irot do 1000 loop
               endif            ! end of if(exchange)
!! ============== End of exchange section ======================
               if(timemix) call timeshow("33333 k-cycle")
!TIME1 "end of exchange section"
!TIME0 


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


               allocate( zwz0(        ntqxx,ntqxx,nstate))
               allocate( zwz (ntqxx,ntqxx,niw*npm,nstate)) 

               zwz0 = 0.0d0
               zwz  = 0.0d0

               allocate( zw (nblochpmx,nblochpmx))
               allocate( zw_zmel(ngb,mpi__Ps:mpi__Pe) )

               ix = 1 + (0 - nw_i) !at omega=0 ! nw_i=0 (Time reversal) or nw_i =-nw
               read(ifrcw,rec=ix) zw ! direct access read Wc(0) = W(0) - v
               zw(:,:) = ( zw(:,:) + dconjg(transpose(zw(:,:))) )/2d0

               do it=mpi__Ts, mpi__Te ! T parallelization
c$$$           do it=1,nstate 
#ifdef USE_GEMM_FOR_SUM
                  call ZGEMM( "C", "N", ngb, mpi__Pn, ngb, dcmplx(1d0,0d0), 
     &                 zw(1,1), nblochpmx, 
     &                 zmel(1,it,mpi__Ps), ngb*nstate,
     &                 dcmplx(0d0,0d0),
     &                 zw_zmel(1,mpi__Ps), ngb )

                  call ZGEMM( "C", "N", ntqxx, mpi__Pn, ngb, dcmplx(1d0,0d0), 
     &                 zmel(1,it,1), ngb*nstate,
     &                 zw_zmel(1,mpi__Ps), ngb,
     &                 dcmplx(0d0,0d0),
     &                 zwz0(1,mpi__Ps,it), ntqxx )
#else
! not USE_GEMM_FOR_SUM
                  do itpp=mpi__Ps, mpi__Pe ! P parallelization
c$$$                 do itpp=1,ntqxx

                     do igb1=1, ngb
                        zw_zmel(igb1,itpp) = sum( dconjg(zw(1:ngb,igb1)) * zmel(1:ngb,it,itpp) )
                     enddo

                     do itp=1, ntqxx
                        zwz0(itp,itpp,it) = sum( dconjg(zmel(1:ngb,it,itp)) * zw_zmel(1:ngb,itpp) )
                     enddo ! itp
                  enddo ! itpp
#endif
               end do ! it

c$$$               write(*,*) "DEBUG: zwz0=",sum(abs(zwz0))

               do ix=1,niw !niw is usually ~10 points.
                  read(ifrcwi,rec=ix) zw ! direct access read Wc(i*omega)=W(i*omega)-v
                  zw(:,:) = ( zw(:,:) + dconjg(transpose(zw(:,:))) )/2d0


                  do it=mpi__Ts, mpi__Te ! T parallelization
c$$$              do it=1,nstate
#ifdef USE_GEMM_FOR_SUM
                     call ZGEMM( "C", "N", ngb, mpi__Pn, ngb, dcmplx(1d0,0d0), 
     &                    zw(1,1), nblochpmx, 
     &                    zmel(1,it,mpi__Ps), ngb*nstate,
     &                    dcmplx(0d0,0d0),
     &                    zw_zmel(1,mpi__Ps), ngb )

                     call ZGEMM( "C", "N", ntqxx, mpi__Pn, ngb, dcmplx(1d0,0d0), 
     &                    zmel(1,it,1), ngb*nstate,
     &                    zw_zmel(1,mpi__Ps), ngb,
     &                    dcmplx(0d0,0d0),
     &                    zwz(1,mpi__Ps,ix,it), ntqxx )
#else
! not USE_GEMM_FOR_SUM
                     do itpp=mpi__Ps, mpi__Pe ! P parallelization
c$$$                 do itpp=1,ntqxx
                        do igb1=1, ngb
                           zw_zmel(igb1,itpp) = sum( dconjg(zw(1:ngb,igb1)) * zmel(1:ngb,it,itpp) )
                        end do

                        do itp=1, ntqxx
                           zwz(itp,itpp,ix,it) = sum( dconjg(zmel(1:ngb,it,itp)) * zw_zmel(1:ngb,itpp) )
                        end do ! itp
                     end do ! itpp
#endif
                  end do ! it


                  if( npm==2 ) then ! w(iw) - w(-iw) Anti Hermitian part
                     read(ifrcwi,rec=ix) zw ! direct access read Wc(i*omega)=W(i*omega)-v
                     zw(:,:) = ( zw(:,:) - dconjg(transpose(zw(:,:))) )/2d0/img

                     do it=mpi__Ts, mpi__Te ! T parallelization
c$$$                 do it=1,nstate
#ifdef USE_GEMM_FOR_SUM
                        call ZGEMM( "C", "N", ngb, mpi__Pn, ngb, dcmplx(1d0,0d0), 
     &                       zw(1,1), nblochpmx, 
     &                       zmel(1,it,mpi__Ps), ngb*nstate,
     &                       dcmplx(0d0,0d0),
     &                       zw_zmel(1,mpi__Ps), ngb )
                        
                        call ZGEMM( "C", "N", ntqxx, mpi__Pn, ngb, dcmplx(1d0,0d0), 
     &                       zmel(1,it,1), ngb*nstate,
     &                       zw_zmel(1,mpi__Ps), ngb,
     &                       dcmplx(0d0,0d0),
     &                       zwz(1,mpi__Ps,ix+niw,it), ntqxx )
#else
! not USE_GEMM_FOR_SUM
                        do itpp=mpi__Ps, mpi__Pe ! P parallelization
c$$$                    do itpp=1,ntqxx
                           do igb1=1, ngb
                              zw_zmel(igb1,itpp) = sum( dconjg(zw(1:ngb,igb1)) * zmel(1:ngb,it,itpp) )
                           enddo
                           
                           do itp=1, ntqxx
                              zwz(itp,itpp,ix+niw,it) = sum( dconjg(zmel(1:ngb,it,itp)) * zw_zmel(1:ngb,itpp) )
                           enddo ! itp
                        enddo ! itpp
#endif
                     end do ! it
                  end if !npm==2
               enddo ! ix

               deallocate( zw )
               deallocate( zw_zmel )

c$$$               write(*,*) "DEBUG: zwz=",sum(abs(zwz))

!! Integration along imag axis for zwz(omega) for given it,itp,itpp
!! itp  : left-hand end of expternal band index.
!! itpp : right-hand end of expternal band index.
!! it   : intermediate state of G.

               allocate( work(ntqxx) )
               do itpp=mpi__Ps, mpi__Pe ! P parallelization
c$$$           do itpp=1,ntqxx

                  work(1:ntqxx) = 0.0d0
                  do itp=1, ntqxx
                     if((jobsw==5).and.(itpp/=itp)) cycle
                     if (jobsw==1.or.jobsw==4) then
                        omegat = ef
                     else
                        omegat = omega(itp)
                     endif

                     do it=mpi__Ts, mpi__Te ! T parallelization
c$$$                 do it=1, nstate
                        we =.5d0*( omegat -ekc(it))
                        if(it <= nctot) then
                           esmrx = 0d0
                        else
                           esmrx = esmr
                        endif

!! ua_auto may be recovered in future...
!! Gaussian smearing. Integration along im axis. zwz(1:niw) and zwz0 are used.
                        work(itp) =  work(itp) +
     &                       wintzsg_npm (npm, zwz(itp,itpp,1:npm*niw,it), zwz0(itp,itpp,it),
     &                       freqx,wx,ua_,expa_,we,niw,esmrx)
                     end do ! it
                  end do ! itp

                  call MPI__AllreduceSumT(work,ntqxx)

                  do itp=1, ntqxx
                     zsec(itp,itpp,ip) = zsec(itp,itpp,ip) + wtt*work(itp)
                  end do ! itp
               end do ! itpp

               deallocate(work)

               deallocate(zwz) !zwzs
               deallocate(zwz0) !zwzs

               if(jobsw==4) goto 2002

!! -------------------------------------------------------------------
!!  Contribution to SEc(qt,w) from the poles of G (integral along real axis)
!!    Currently, jobsw =1,3,5 are allowed...
!!    The variable we means \omega_epsilon in Eq.(55) in PRB76,165106 (2007)
!! -------------------------------------------------------------------
               if(timemix) call timeshow("goto Sec pole part k-cycle")
               if(debug)  write(6,*)'GOTO contribution to SEc(qt,w) from the poles of G'
               if (.not.(jobsw == 1 .or. jobsw == 3.or.jobsw==5)) then
                  call rx( 'sxcf_fal3_scz: jobsw /= 1 3 5')
               endif
!! Get index nwxi nwx nt_max. finish quickly. We can simplify this...
               call get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,
     i              nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,
     o              nwxi,nwx,nt_max)
!! assemble small arrays first.
               allocate(we_(nt_max,ntqxx),wfac_(nt_max,ntqxx),ixss(nt_max,ntqxx),ititpskip(nt_max,ntqxx),iirx(ntqxx))
               call weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,
     i              ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt,
     o              we_,wfac_,ixss,ititpskip,iirx)

!! We need zw3, the Hermitian part, because we need only hermitean part of Sigma_nn'
!! This can be large array; nwx-nwxi+1 \sim 400 or so...
               allocate( zw(nblochpmx,nblochpmx))
               allocate( zw3(ngb,ngb,nwxi:nwx)) 

               do ix = nwxi,nwx
                  nrec = ix-nw_i+1 !freq_r(ix is in nw_i:nx)
                  read(ifrcw,rec=nrec) zw ! direct access Wc(omega) = W(omega) - v
                  zw3(1:ngb,1:ngb,ix)=(zw(1:ngb,1:ngb)+transpose(dconjg(zw(1:ngb,1:ngb))))/2d0
               enddo
               deallocate(zw)

!! jobsw==3
               if( jobsw==3 ) then
                  allocate(zw_zmel(ngb,3))
                  allocate(zmel_zw_zmel(3))
                  allocate(work(ntqxx))

                  do itp=1,ntqxx
                     work(1:ntqxx) = 0.0d0

                     do it=1,nt_max ! T parallelization
                        if( it<mpi__Ts .or. mpi__Te<it ) cycle

                        if(ititpskip(it,itp)) cycle
                        we =  we_(it,itp) 
                        ixs=  ixss(it,itp)

                        do ix0=1,3
                           ix=ixs+ix0-2
#ifdef USE_GEMM_FOR_SUM
                           call ZGEMV( "C", ngb, ngb, dcmplx(1d0,0d0), 
     &                          zw3(1,1, iirx(itp)*ix), ngb,
     &                          zmel(1,it,itp), 1,
     &                          dcmplx(0d0,0d0),
     &                          zw_zmel(1,ix0), 1 )
#else
! not USE_GEMM_FOR_SUM
                           do igb2=1,ngb 
                              zw_zmel(igb2,ix0) = sum( dconjg(zw3(1:ngb,igb2, iirx(itp)*ix)) * zmel(1:ngb,it,itp) )
                           end do ! igb2
#endif
                        end do ! ix0

                        do itpp=mpi__Ps, mpi__Pe ! P parallelization
c$$$                    do itpp=1,ntqxx
                           do ix0=1,3
#ifdef USE_GEMM_FOR_SUM
                              zmel_zw_zmel(ix0) = ZDOTC( ngb, zw_zmel(1,ix0), 1, zmel(1,it,itpp), 1 )
#else
! not USE_GEMM_FOR_SUM
                              zmel_zw_zmel(ix0) = sum( dconjg(zw_zmel(1:ngb,ix0)) * zmel(1:ngb,it,itpp) )
#endif
                           end do ! ix0

c$$$                           if(npm==1) then
c$$$                              zsec(itp,itpp,ip) = zsec(itp,itpp,ip) + wfac_(it,itp)
c$$$     &                             * alagr3z2(we,freq_r(ixs-1),zmel_zw_zmel,itp==itpp )
c$$$                           else
c$$$                              zsec(itp,itpp,ip) = zsec(itp,itpp,ip) + wfac_(it,itp)
c$$$     &                             * alagr3z(we,freq_r(ixs-1),zmel_zw_zmel)
c$$$                           endif

                           if(npm==1) then
                              work(itpp) = work(itpp) + wfac_(it,itp)
     &                             * alagr3z2(we,freq_r(ixs-1),zmel_zw_zmel,itp==itpp )
                           else
                              work(itpp) = work(itpp) + wfac_(it,itp)
     &                             * alagr3z(we,freq_r(ixs-1),zmel_zw_zmel)
                           endif
                        end do ! itpp
                     end do ! it

                     call MPI__AllreduceSumT(work,ntqxx)
                     
                     do itpp=mpi__Ps, mpi__Pe ! P parallelization
c$$$                 do itpp=1,ntqxx
                        zsec(itp,itpp,ip) = zsec(itp,itpp,ip) + work(itpp)
                     enddo ! itpp
                  enddo ! itp

                  deallocate(work)
                  deallocate(zmel_zw_zmel)
                  deallocate(zw_zmel)

               endif ! jobsw==3

!! jobsw=1,5 Sigma are calculated.
               if( jobsw==1.or.jobsw==5) then
                  allocate(zw_zmel(ngb,3))
                  allocate(zmel_zw_zmel(3))
                  allocate(work(ntqxx))

                  do itp=mpi__Ps, mpi__Pe ! P parallelization
c$$$              do itp=1,ntqxx
                     work(1:ntqxx) = 0.0d0

                     do it=1,nt_max ! T parallelization
                        if( it<mpi__Ts .or. mpi__Te<it ) cycle

                        if(ititpskip(it,itp)) cycle
                        we =  we_(it,itp) 
                        ixs=  ixss(it,itp)
                        zwz3=0d0

                        do ix0=1,3
                           ix=ixs+ix0-2
#ifdef USE_GEMM_FOR_SUM
                           call ZGEMV( "N", ngb, ngb, dcmplx(1d0,0d0), 
     &                          zw3(1,1,iirx(itp)*ix), ngb,
     &                          zmel(1,it,itpp), 1,
     &                          dcmplx(0d0,0d0),
     &                          zw_zmel(1,ix0), 1 )
                           zmel_zw_zmel(ix0) = ZDOTC( ngb, zmel(1,it,itp), 1, zw_zmel(1,ix0), 1 )
#else
! not USE_GEMM_FOR_SUM
                           do igb1=1, ngb
                              zw_zmel(igb1,ix0) = sum( zw3(igb1,1:ngb,iirx(itp)*ix) * zmel(1:ngb,it,itp) )
                           end do
                           zmel_zw_zmel(ix0) = sum( dconjg(zmel(1:ngb,it,itp)) * zw_zmel(1:ngb,ix0) )
#endif

c$$$!!               **** most time consuming part for jobsw=1 ******
c$$$!!               To reduce computational time, confusing treatment only uses lower half of zw3 (zw3 is Hermitan)
c$$$!!               Clean up needed.
c$$$                           do igb2=2,ngb 
c$$$                              zz2 = sum( dconjg(zmel(1:igb2-1,it,itp)) *zw3(1:igb2-1,igb2,iirx(itp)*ix)  ) +
c$$$     &                             .5d0* dconjg(zmel(igb2,it,itp)) *zw3(igb2,igb2,iirx(itp)*ix)
c$$$                              zwz3(ix0) = zwz3(ix0)+zz2*zmel(igb2,it,itp)
c$$$                           end do !igb2
c$$$                           zwz3(ix0)=2d0*dreal(zwz3(ix0))+ !I think 2d0 is from upper half.
c$$$     &                          dconjg(zmel(1,it,itp))*zw3(1,1, iirx(itp)*ix)*zmel(1,it,itp)



!     a little complicated but zwz3 contains     <itp| it I> wz3_IJ(we)  <J it| itp>
                        end do ! ix0

c$$$                        if(npm==1) then
c$$$                           zsec(itp,itp,ip) = zsec(itp,itp,ip) 
c$$$     &                          + wfac_(it,itp)*alagr3z2(we,freq_r(ixs-1),zmel_zw_zmel,itp==itp)
c$$$                        else
c$$$                           zsec(itp,itp,ip) = zsec(itp,itp,ip) 
c$$$     &                          + wfac_(it,itp)*alagr3z(we,freq_r(ixs-1),zmel_zw_zmel)
c$$$                        endif

                        if(npm==1) then
                           work(itp) = work(itp) 
     &                          + wfac_(it,itp)*alagr3z2(we,freq_r(ixs-1),zmel_zw_zmel,itp==itp)
                        else
                           work(itp) = work(itp) 
     &                          + wfac_(it,itp)*alagr3z(we,freq_r(ixs-1),zmel_zw_zmel)
                        endif
                     end do ! it
                  enddo ! itp

                  call MPI__AllreduceSumT(work,ntqxx)

                  do itp=mpi__Ps, mpi__Pe ! P parallelization
c$$$              do itp=1,ntqxx
                     zsec(itp,itp,ip) = zsec(itp,itp,ip) + work(itp)
                  end do ! itp

               end if            ! this contribution to zsec_nn is real (hermitean)

               deallocate( zw3 )
               goto 2012

 2012          continue
               deallocate(we_,wfac_,ixss,ititpskip,iirx)
 2002          continue
               if(verbose()>50) call timeshow("11after alagr3z iw,itp,it cycles")
               if(debug) then
                  write(6,*)' end of do 2001 '
                  do itp=mpi__Ps, mpi__Pe ! P parallelization
c$$$              do itp=1,ntq
                     write(6,'(" zsec=",i3,2d15.7)') itp,zsec(itp,itp,ip)
                  enddo
               endif
!TIME1 "end of one of do doop 1000 iteration"
!TIME0 
 1000       continue            ! end do irot
            ifvcoud =iclose('Vcoud.'//charnum5(kx))
            if(.not.exchange) then
               ifrcw  = iclose('WVR.'//charnum5(kx))
               ifrcwi = iclose('WVI.'//charnum5(kx))
            endif   
!TIME1 "end of 1100 loop"
            deallocate(ppovlz)
 1100    continue               ! end of kx-loop

         ifvcoud =iclose('Vcoud')
         if(irot==1) write(6,"('  sum(abs(zsec))=',d23.15)") sum(abs(zsec))
         if (allocated(vcoul))deallocate(vcoul)
 1001 continue                  ! end do ip
!TIME1 "end of sxcf_fal2_scz"


      end subroutine sxcf_fal3_scz



      subroutine weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,
     i     ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt,
     o     we_,wfac_,ixss,ititpskip,iirx)
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
!     jobsw==2  
!     if (jobsw==2)  omg=.5d0*(omega(itp)+omega(itpp))
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
         do it = itini,itend    ! nt0p corresponds to efp
            esmrx = esmr
            if(it<=nctot) esmrx = 0d0
            wfac_(it,itp) = wfacx2(omg,ef, ekc(it),esmrx)
            wfac = wfac_(it,itp)
            if(wfac<wfaccut) then
               ititpskip(it,itp)=.true.
               cycle 
            endif  
            wfac_(it,itp)=  wfac_(it,itp)*wtt*iii
!     Gaussian smearing we_= \bar{\omega_\epsilon} in sentences next to Eq.58 in PRB76,165106 (2007)
!     wfac_ = $w$ weight (smeared thus truncated by ef). See the sentences.
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

      end module m_sxcf



!! --------------------------------------------------------------
      subroutine get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,
     i     nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,
     o     nwxi,nwx,nt_max)
!     > Determine indexes of a range for calculation.
!! It is better to clean this up...
      implicit none
      integer,intent(in) :: nctot,nw_i,nw,nstate,nt0p,nt0m,ntq,
     &     nband,ntqxx
      real(8),intent(in):: omega(ntq),esmr,ef,ekc(nctot+nband),wfaccut,
     &     freq_r(nw_i:nw)
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
            if(it<=nctot) then
               if(wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
            endif
            do iwp = 1,nw
               ixs=iwp
               if(freq_r(iwp)>abs(we)) exit
            enddo
            if(ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
            if(ixs>ixsmin .and. omg< ef ) ixsmin = ixs
            wexx  = we
            if(ixs+1 > nw) then
               write (*,*) ' nw_i ixsmin',nw_i, ixsmin
               write (*,*) ' wexx ',wexx
               write (*,*) ' omg ekc(it) ef ', omg,ekc(it),ef
               call rx( ' sxcf 222: |w-e| out of range')
            endif
 311     continue
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
               if(wfac>wfaccut) then
                  if (it > nt_max) nt_max=it ! nt_max is  unocc. state
               endif            ! that ekc(it>nt_max)-omega > 0
            enddo               ! so it > nt_max does not contribute to omega pole integral
         endif
 401  continue                  !end of  w and qt -loop
      end subroutine get_nwx
