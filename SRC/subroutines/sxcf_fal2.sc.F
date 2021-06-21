!> this module is only because name=name argument binding. No data
      module m_sxcfsc
      use m_readqg,only   : readqg0
      use m_readeigen,only: readeval
      use m_keyvalue,only   : getkeyvalue
      use m_zmel,only     : Get_zmel_init, Get_zmel,Get_zmeltt,Setppovlz, Setppovlz_chipm, Deallocate_zmel,
     o     zmel,zmeltt
      use m_itq,only: itq,ntq
      use m_genallcf_v3,only: nlmto,nsp=>nspin,nctot,ngrp,niw,ecore !,symgg
      use m_READ_BZDATA,only: qbas,ginv,nstbz,qibz,qbz,wk=>wbz,nqibz,nqbz
c      use m_freq,only: freq_r,nw_i,nw, freqx,wx !,dwdummy,         
      use m_readVcoud,only:   Readvcoud, vcoud,vcousq,zcousq,ngb,ngc
      use m_readfreq_r,only: freq_r, nw_i,nw,freqx,wx=>wwx,nblochpmx
      use m_rdpp,only: Rdpp,    !"call rdpp" generate following data.
     &     nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
      use m_readqg,only:  ngpmx,ngcmx
      use m_readq0p,only: ReadQ0P, nq0i, wqt,q0i
      use m_readhbe,only: nband,mrecg
      use m_readepswklm,only: wklm,lxklm
      use m_readfreq_r,only: mrecl,expa_,npm
      implicit none
!---------------------------------      
      public sxcf_scz
      contains
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      subroutine sxcf_scz(qip,ef,esmr,isp,irkip,nrkip,nq,exchange,jobsw,hermitianW,
     o zsec,coh,nbandmx)
      intent(in)          qip,ef,esmr, 
     i isp, irkip,nrkip,  nq,  !wgt0,
     i exchange,!screen, 
     i jobsw, hermitianW
     i          ,nbandmx
      intent(out) :: zsec , coh
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
!!
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
!!  1,3,5scGW mode.
!!   diag+@EF      jobsw==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
!!   xxx modeB (Not Available now)  jobsw==2 SE_nn'((e_n+e_n')/2)  !we need to recover comment out for jobsw==2, and test.
!!   mode A        jobsw==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2 (Usually usued in QSGW).
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
!! -------------------------------------------------------------------------------
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
!! -------------------------------------------------------------------------------
!! \endverbatim
!! \verbatim
!!
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
      logical :: exchange!,screen!,cohtest
      integer :: nq !,natom ntq,nqbz,nqibz,
      integer :: isp !nband,nq0i,isp!,nsp !,mdim(*) !,nlnmx nlmto,nctot,
      integer :: irkip(nqibz,ngrp,nq),nrkip(nqibz,ngrp,nq)
      integer :: jobsw !kount(nqibz,nq),ngpmx,ngcmx,ifexsp,
c      integer :: lxklm !,invg(ngrp) !nnr(nkpo),nor(nkpo),
      real(8) :: ef,esmr!,wklm((lxklm+1)**2) !,qrr(3,nkpo),dwdummy freqx(niw),wx(niw),
      real(8) :: qip(3,nq)!,eftrue!,ecorem !,ebmx
      integer :: nbandmx(nq)
      complex(8),optional :: zsec(ntq,ntq,nq) , coh(ntq,nq)
      real(8):: ebmx
      integer :: ifrcw,ifrcwi
      integer :: ip, it, itp, i, ix, kx, irot, kr
      integer :: nt0p, nt0m,nstate , nbmax, ntqxx !iatomp(natom), 
      integer :: nt,ixs,iw,ivc,ifvcoud,ngb0
      integer :: nprecx,mrecl,ifwd,nrot,nwp,ierr 
      integer :: nstatetot,iqini,iqend!, ngb,ngc !nbcut,
      integer :: invr,ia,nn,ntp0,no,itpp,nrec,npm,itini,itend,nbmxe
      integer :: iwp,nwxi,nwx,iir, igb1,igb2,ix0,iii

      real(8) :: tpi, ekc(nctot+nband),ekq(nband), det, q(3),ua_
      real(8) :: expa_(niw), qxx(3), symope(3,3),shtv(3) !tr(3,natom), 
      real(8) :: efp,efm,wtt,wfac,we,esmrx,qbasinv(3,3)
      real(8) :: qvv(3),pi,fpi,eq(nband),omega(ntq),quu(3),freqw,ratio
      real(8) :: qibz_k(3),qbz_kr(3),ddw,vc,omega0,omg

c      complex(8) :: cphiq(nlmto,nband), cphikq(nlmto,nband)
      complex(8) :: zwzs0,zz2,zwz3(3)

! local arrays
c      real(8),intent(in) :: freq_r(nw_i:nw)
      real(8),allocatable :: drealzzzmel(:,:,:), dimagzzzmel(:,:,:),uaa(:,:)
      complex(8),allocatable :: vcoul(:,:) 
      complex(8),allocatable :: zzzmel(:,:,:)
      complex(8),allocatable :: zw (:,:)
      complex(8),allocatable:: zwz(:,:,:,:)
      complex(8),allocatable :: zwz0(:,:,:)!,zwzi(:,:,:)
      complex(8),allocatable :: zwix(:,:),zwzix(:,:,:),zmel1(:) !,expikt(:)
      complex(8), allocatable :: zmel1_(:,:,:), zw3(:,:,:),zw3x(:,:)
      complex(8), allocatable :: zwz4(:,:),zwz44(:,:),pomat(:,:), zwzs(:)
      complex(8),allocatable :: ppovl(:,:),zcousq(:,:)

      real(8), parameter :: wfaccut=1d-8
      complex(8), parameter :: img=(0d0,1d0)

! external function
c      logical :: smbasis
c      logical :: test_symmetric_W
c      logical :: GaussSmear !fixed to be T
c      logical :: newaniso !fixed to be T
c      integer :: bzcase !fixed to be 1
      character(5) :: charnum5
c      integer :: iopen,iclose
      integer :: invrot
c      complex(8) ::  wintzsg_npm !wintzav,
      integer :: nocc
      real(8) :: wfacx
      real(8) :: wfacx2
      real(8) :: weavx2
      complex(8) :: alagr3z
      complex(8) :: alagr3z2,zadd

      integer:: ndummy1,ndummy2,nlmtobnd,nt0
      real(8):: wexx
c      complex(8),allocatable :: z1p(:,:,:),vcoult(:,:)
      logical :: debug, debugp,debug2=.false.
c      logical :: gass           !external
c      real(8):: wgtq0p
      integer::verbose,ififr !,ifile_handle
      real(8):: ua2_(niw),freqw1
      integer :: istate,  nt_max !nbcutc,nbcutin, 
      real(8):: q_r(3),qk(3)!,omegat
      logical::  oncew, onceww, eibz4sig,  timemix

      integer,allocatable:: ixss(:,:),iirx(:)
      real(8),allocatable:: we_(:,:),wfac_(:,:)
      complex(8),allocatable:: zw3av(:,:),zmelw(:,:,:),w3p(:),wtff(:)
      integer:: noccx
      real(8)::polinta
      logical,allocatable:: ititpskip(:,:)

      logical:: tote=.false.
      logical:: hermitianW

      real(8),allocatable:: wcorehole(:,:),vcoud_(:),wfft(:)
      logical:: corehole
      integer:: ifcorehole,ixx
      real(8):: tolq=1d-8!,wgt3(3)
!!----------------------------------------------------------------
      if (.not.(jobsw==1.or.jobsw==3.or.jobsw==5)) call rx('sxcf_scz: jobsw=1or3or5')
      if(verbose()>=90) debug=.true.
      timemix=.false.
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      debug=.false.
!! core-hole mode can work but we need to turn it on by hand.
      corehole=.false.
      if(corehole) then
        open(newunit=ifcorehole,file='CoreHole') 
        if(allocated(wcorehole)) deallocate(wcorehole)
        allocate(wcorehole(nctot,nsp))
        do it=1,nctot
          read(ifcorehole,*) wcorehole(it,1:nsp)
        enddo
        close(ifcorehole)
        write(*,*) 'end of reading CoreHole' 
      endif
      if(.not.exchange) then
        open(newunit=ifwd,file='WV.d')
        read(ifwd,*) nprecx,mrecl
        close(ifwd) 
!! gauss_img : interpolation gaussion for W(i \omega).
        call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
        do ix = 1,niw           !! Energy mesh; along im axis.
          freqw     = (1d0 - freqx(ix))/ freqx(ix)
          expa_(ix) = exp(-(ua_*freqw)**2)
        enddo
        npm = 1                 ! npm=1    Timeveversal case
        if(nw_i/=0) npm = 2     ! npm=2 No TimeReversal case. Need negative energy part of W(omega)
      endif 
      tpi         = 8d0*datan(1d0)
      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot,isp)!-ecorem ! core
      nlmtobnd    = nlmto*nband
      nstatetot   = nctot + nband

!!== ip loop to spedify external q ==
      do 1001 ip = 1,nq         
        if(sum(irkip(:,:,ip))==0) cycle ! next ip
        write (6,*) ip,'  out of ',nq,'  k-points(extrnal q) '
        q(1:3)= qip(1:3,ip)
        eq= readeval(q,isp)
        do i  = 1,ntq
          omega(i) = eq(itq(i)) 
        enddo
!! we only consider bzcase()==1
        if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
        if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
!! NOTE total number of 
!!    kx loop(do 1100) and irot loop (do 1000) makes all the k mesh points.
!!    When iqini=1 (Gamma point), we use effective W(q=0) defined in the paper.
        iqini = 1
        iqend = nqibz             !no sum for offset-Gamma points.
        do 1100 kx = iqini,iqend 
          if(sum(irkip(kx,:,ip))==0) cycle ! next kx
          write(6,*) ' ### do 1100 start kx=',kx,' from ',iqini,' through', iqend
          qibz_k = qibz(:,kx)
          if(timemix) call timeshow("11111 k-cycle")
          call Readvcoud(qibz_k,kx,NoVcou=.false.) ! Readin ngc,ngb,vcousq,zcousq !for the Coulomb matrix
          call Setppovlz(qibz_k,matz=.true.)
          if(debug) write(6,*) ' sxcf_fal2sc: ngb ngc nbloch=',ngb,ngc,nbloch
          
!! === open WVR,WVI ===
          if(.not.exchange) then
            open(newunit=ifrcw,  file='WVR.'//charnum5(kx), action='read',form='unformatted',
     &          status='unknown',access='direct',recl=mrecl)
            open(newunit=ifrcwi, file='WVI.'//charnum5(kx), action='read',form='unformatted',
     &          status='unknown',access='direct',recl=mrecl)
          endif
          nrot=0
          do irot = 1,ngrp
             kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
             if(kr==0) cycle    ! next irot
             qbz_kr= qbz (:,kr) 
             nrot=nrot+1
          enddo  

!! === loop 1000 over rotations irot ===
          do 1000 irot = 1,ngrp
            kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
            if(kr==0) cycle
            qbz_kr= qbz (:,kr) 
!! no. occupied (core+valence) and unoccupied states at q-rk
            qk =  q - qbz_kr        
            ekq = readeval(qk, isp) 
            ekc(nctot+1:nctot+nband) = ekq (1:nband)
            nt0 = nocc (ekc,ef,.true.,nstatetot)
            ddw= 10d0
            efp= ef+ddw*esmr
            efm= ef-ddw*esmr
            nt0p = nocc (ekq,efp,.true.,nstatetot)+ nctot
            nt0m = nocc (ekq,efm,.true.,nstatetot)+ nctot
            if(exchange) then
              nbmax = nt0p-nctot
            else
              ebmx=1d10 !this is needed probably because to fill 1d99 for ekc(i) above boundary.
              nbmxe = nocc (ekc,ebmx,.true.,nstatetot)-nctot
              nbmax  = min(nband,nbmxe) !nbmx_sig,
            endif
!! ntqxx is number of bands for <i|sigma|j>.
            ntqxx = nbandmx(ip) !mar2015
            if(debug) write(6,*)' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
            nstate = nctot + nbmax ! = nstate for the case of correlation
!! Get matrix element zmelt= rmelt + img*cmelt, defined in m_zmel.F---
! this return zmeltt (for exchange), or zmel (for correlation)
            call Get_zmel_init(exchange,q,qibz_k,irot,qbz_kr,isp, nbmax,ntqxx,nctot,ncc=0)
            if(     exchange) call Get_zmeltt()!return zmeltt 
            if(.not.exchange) call Get_zmel()  !rerutn zmel
            if(debug) write(6,*)' end of get_zmelt'
            wtt = wk(kr)
            if(eibz4sig()) then
              wtt=wtt*nrkip(kx,irot,ip) 
            endif   
!!--------------------------------------------------------
!! --- exchange section ---
!!--------------------------------------------------------
            if(exchange) then   !At the bottom of this block, cycle do 1000 irot.
c We need to check following comments ----     
c     S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c     v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
c     
c     > z1p(j,n,t) = S[i=1,nbloch] <psi(q,t) | psi(q-rk,n) B(rk,i)> v(k)(i,j)
              allocate(vcoud_(ngb),wtff(nctot+nbmax),w3p(nctot+nbmax)) 
              vcoud_= vcoud
              if(kx == iqini) vcoud_(1) = wklm(1)* fpi*sqrt(fpi) /wk(kx)
              do it=1,nctot+nbmax
                 wtff(it)= wfacx(-1d99, ef, ekc(it), esmr)
              enddo   
              do 8010 itpp= 1,ntqxx
              do 8011 itp = 1,ntqxx !S[j=1,nbloch]  z1p(j,t,n) <B(rk,j) psi(q-rk,n) |psi(q,t')>
                if(jobsw==5.and.(itpp/=itp)) cycle
                do it = 1, nctot+nbmax
                   w3p(it)= sum(zmeltt(it,itp,:) * vcoud_(:)* dconjg(zmeltt(it,itpp,:)))
                enddo
                w3p(nctot+1:nctot+nbmax) = w3p(nctot+1:nctot+nbmax) * wtff(nctot+1:nctot+nbmax)
                if(corehole) w3p(1:nctot) = w3p(1:nctot) * wcorehole(1:nctot,isp) 
                zsec(itp,itpp,ip) = zsec(itp,itpp,ip) - wtt * sum( w3p(:) )
 8011         continue
 8010         continue
              deallocate(vcoud_,wtff,w3p)
              cycle             ! next irot do 1000 loop
            endif               ! end of if(exchange)
!! ============== End of exchange section ======================
            if(timemix) call timeshow("33333 k-cycle")
            
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
c            allocate( zwz0(        nstate,ntqxx,ntqxx))
            allocate( zwz(0:niw*npm,nstate,ntqxx,ntqxx)) 
            allocate( zw (nblochpmx,nblochpmx))
            do ixx=0,niw        !niw is usually ~10 points.
               if(ixx==0) then
                  ix = 1 + (0 - nw_i) !at omega=0 ! nw_i=0 (Time reversal) or nw_i =-nw
                  read(ifrcw,rec=ix) zw ! direct access read Wc(0) = W(0) - v
               elseif(ixx>0) then
                  ix = ixx
                  read(ifrcwi,rec=ix) zw ! direct access read Wc(i*omega)=W(i*omega)-v
               endif
               call matzwz2(2, zw(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb, zwz(ixx,1:nstate,1:ntqxx,1:ntqxx)) ! zwz = zmel*(W(i w)-v)*zmel
               do istate = 1,nstate
                  block
                  complex(8):: zwx(ntqxx,ntqxx)
                  zwx = zwz(ix,istate,1:ntqxx,1:ntqxx)
                  zwz(ix,istate,1:ntqxx,1:ntqxx) = ( zwx+dconjg(transpose(zwx)))/2d0                   ! w(iw) + w(-iw)  Harmitian part
                  if(npm==2) zwz(ix+niw,istate,1:ntqxx,1:ntqxx) = (zwx-dconjg(transpose(zwx)))/2d0/img ! w(iw) - w(-iw) Anti Hermitian part
                  end block
               enddo   
            enddo
            deallocate(zw)
!! Integration along imag axis for zwz(omega) for given it,itp,itpp
!! itp  : left-hand end of expternal band index.
!! itpp : right-hand end of expternal band index.
!! it   : intermediate state of G.
c            allocate(zwzi(ntqxx,ntqxx))
            block
            real(8):: esmrx(nstate),omegat(ntqxx)
            real(8):: wgtim(0:npm*niw,ntqxx,nstate)
!     , wintzsg_npm_wgtim(0:npm*niw),wintzsg_npm_wgtimxxx(0:npm*niw)
            esmrx(1:nctot)=0d0
            esmrx(nctot+1:nstate)=esmr
            omegat(1:ntqxx) = omega(1:ntqxx)
            if (jobsw==1) omegat = ef
!     ! Gaussian smearing. Integration weight wgtim along im axis for zwz(0:niw*npm) 
            do itp = 1,ntqxx
            do  it = 1,nstate
               we = .5d0*(omegat(itp)-ekc(it))
               call wintzsg_npm_wgtim(npm,freqx,wx,ua_,expa_,we,niw,esmrx(it), wgtim(:,itp,it))
            enddo
            enddo
            do 1400 itpp= 1,ntqxx
            do 1410 itp = 1,ntqxx
               if((jobsw==5).and.(itpp/=itp)) cycle
               zsec(itp,itpp,ip) = zsec(itp,itpp,ip) + wtt*sum(wgtim(:,itp,:)*zwz(:,:,itp,itpp))
 1410       continue  
 1400       continue
            end block
            deallocate(zwz) 

!! -------------------------------------------------------------------
!!  Contribution to SEc(qt,w) from the poles of G (integral along real axis)
!!    Currently, jobsw =1,3,5 are allowed...
!!    The variable we means \omega_epsilon in Eq.(55) in PRB76,165106 (2007)
!! -------------------------------------------------------------------
            if(debug)  write(6,*)'GOTO contribution to SEc(qt,w) from the poles of G'
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
c     allocate(zwz4(ntqxx,3))
              block
              real(8)::wgt3(3,nt_max,ntqxx)
              do itp=1,ntqxx
              do it=1,nt_max
                  if(ititpskip(it,itp)) cycle
                  ixs=  ixss(it,itp)
                  we =  we_(it,itp) 
                  if(npm==1) call alagr3z2wgt(we,freq_r(ixs-1),wgt3(:,it,itp))
                  if(npm==2) call alagr3zwgt (we,freq_r(ixs-1),wgt3(:,it,itp))
               enddo
               enddo
               do itp=1,ntqxx
               do it=1,nt_max
                  if(ititpskip(it,itp)) cycle
                  ixs=  ixss(it,itp)
                  do itpp=1,ntqxx
                     zadd =    wgt3(1,it,itp)* sum(matmul( dconjg(zmel(:,it,itp)),
     &                          zw3(:,:,iirx(itp)*(ixs+1-2)))*zmel(:,it,itpp))
     &                       + wgt3(2,it,itp)* sum(matmul( dconjg(zmel(:,it,itp)),
     &                          zw3(:,:,iirx(itp)*(ixs+2-2)))*zmel(:,it,itpp))
     &                       + wgt3(3,it,itp)* sum(matmul( dconjg(zmel(:,it,itp)),
     &                          zw3(:,:,iirx(itp)*(ixs+3-2)))*zmel(:,it,itpp))
                     if(npm==1) then
                        if(itp==itpp.and.dimag(zadd)>0d0) zadd=dreal(zadd)
                        zsec(itp,itpp,ip)=zsec(itp,itpp,ip)+ wfac_(it,itp)* zadd
                     else
                        zsec(itp,itpp,ip)=zsec(itp,itpp,ip)+ wfac_(it,itp)* zadd
                     endif
                  enddo
                enddo
              enddo
              end block
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
!!                  **** most time consuming part for jobsw=1 ******
!!                  To reduce computational time, confusing treatment only uses lower half of zw3 (zw3 is Hermitan)
!!                  Clean up needed.
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
 2012       continue
            deallocate(we_,wfac_,ixss,ititpskip,iirx)
 2002       continue
            deallocate(zw3, zmel1)
            call Deallocate_zmel()
            if(allocated(zmel1_)) deallocate(zmel1_)
            if(verbose()>50) call timeshow("11after alagr3z iw,itp,it cycles")
            if(debug) then
              write(6,*)' end of do 2001 '
              do itp = 1,ntq
                write(6,'(" zsec=",i3,2d15.7)') itp,zsec(itp,itp,ip)
              enddo
            endif
 1000     continue              ! end do irot
          close(ifvcoud)
          if(.not.exchange) then
            close(ifrcw)
            close(ifrcwi)
          endif   
 1100   continue                ! end of kx-loop
        if(irot==1) write(6,"('  sum(abs(zsec))=',d23.15)") sum(abs(zsec))
        if (allocated(vcoul)) deallocate(vcoul)
 1001 continue                  ! end do ip
      end subroutine sxcf_scz
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      subroutine weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,
     i     ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt,
     o     we_,wfac_,ixss,ititpskip,iirx)
      implicit none
      intent(in)::                 nctot,esmr,omega,ekc,freq_r,nw_i,nw,
     i     ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt
      intent(out)::
     o     we_,wfac_,ixss,ititpskip,iirx
!! generate required data set for main part of real part integration.
      integer:: ntqxx,nctot,nw_i,nw,nt0m,nwx,nwxi,nt_max
      real(8):: ef,omega(ntqxx),ekc(ntqxx),freq_r(nw_i:nw),esmr,wfaccut,wtt
      real(8):: we_(nt_max,ntqxx),wfac_(nt_max,ntqxx)
      integer:: ixss(nt_max,ntqxx),iirx(ntqxx)
      logical:: ititpskip(nt_max,ntqxx)
      integer:: itini,iii,it,itend,wp,ixs,itp,iwp,nt0p
      real(8):: omg,esmrx,wfacx2,we,wfac,weavx2
      ititpskip=.false.
      do itp = 1,ntqxx          !this loop should finish in a second
        omg = omega(itp)
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
      
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
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
