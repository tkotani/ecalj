      module m_sxcfsc !self-energy calculation 
      use m_readqg,only   : Readqg0
      use m_readeigen,only: Readeval
      use m_keyvalue,only   : Getkeyvalue
      use m_zmel,only : Get_zmel_init, Setppovlz, Setppovlz_chipm, Deallocate_zmel, zmel
      use m_itq,only: itq,ntq
      use m_genallcf_v3,only: nlmto,nspin,nctot,niw,ecore !,symgg
      use m_READ_BZDATA,only: qibz,qbz,wk=>wbz,nqibz,nqbz,wklm,lxklm,nq0i, wqt=>wt,q0i
      use m_readVcoud,only:   Readvcoud, vcoud,ngb,ngc
      use m_readfreq_r,only: freq_r, nw_i,nw,freqx,wx=>wwx,nblochpmx,mrecl,expa_,npm,nprecx
      use m_rdpp,only: Rdpp, nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
      use m_readqg,only:  ngpmx,ngcmx
      use m_readhbe,only: nband,mrecg
      use m_hamindex,only: ngrp
      use m_eibzhs,only: nrkip=>nrkip_all,irkip_all
      use m_readgwinput,only: ua_,  corehole,wcorehole
      use m_mpi,only: MPI__sxcf_rankdivider
      implicit none
!---------------------------------      
      public sxcf_scz, zsecall
      complex(8),allocatable,target:: zsecall(:,:,:,:) !output
      private
      contains
      
      subroutine sxcf_scz(qip,ef,esmr,nq,exchange,jobsw,hermitianW,nbandmx,ixc,nspinmx)
      intent(in)          qip,ef,esmr,nq,exchange,jobsw,hermitianW,nbandmx,ixc,nspinmx
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
!!
!! \remark
!! jobsw switch. We now support only mode=3-----------------
!!  1,3,5scGW mode.
!!   diag+@EF      jobsw==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
!!   modeB (Not Available now)  jobsw==2 SE_nn'((e_n+e_n')/2) 
!!   mode A        jobsw==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2 (Usually usued in QSGW).
!!   diagonly      jobsw==5 delta_nn' SE_nn(e_n) (not efficient memoryuse; but we don't use this mode so often).
!!
!! \verbatim
!!  eftrue is added. !Jan2013
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
!!    qlat    = base reciprocal lattice vectors
!!    ginv    = inverse of qbas =transpose(plat)
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
      logical :: exchange
      integer :: nq,isp,nspinmx,jobsw 
      integer :: nbandmx(nq,nspinmx)
      real(8) :: ef,esmr, qip(3,nq)
      real(8):: ebmx
c     complex(8),optional :: zsec(ntq,ntq,nq)
      complex(8),pointer::zsec(:,:)
      
      integer,allocatable :: ifrcw(:),ifrcwi(:)
      integer :: ip, it, itp, i, ix, kx, irot, kr
      integer :: nt0p, nt0m,nstate , nbmax, ntqxx 
      integer :: nt,ixs,iw,ivc,ifvcoud,ngb0
      integer :: ifwd,nrot,nwp,ierr 
      integer :: iqini,iqend
      integer :: invr,ia,nn,ntp0,no,itpp,nrec,itini,itend,nbmxe
      integer :: iwp,nwxi,nwx,iir, igb1,igb2,ix0,iii
      integer :: invrot,nocc,nlmtobnd,nt0,verbose,ififr, istate,  nt_max ,noccx
      real(8) :: ekc(nctot+nband),ekq(nband), det, q(3) !,ua_
      real(8) :: wtt,wfac,we,esmrx
      real(8) :: qvv(3),eq(nband),omega(ntq),quu(3),freqw,ratio
      real(8) :: qibz_k(3),qbz_kr(3),vc,omega0,omg
      complex(8),allocatable:: zwz(:,:,:), zw (:,:)
      real(8), parameter :: wfaccut=1d-8,tolq=1d-8
      complex(8), parameter :: img=(0d0,1d0)
      character(5) :: charnum5
      real(8)::polinta, wfacx, wfacx2, weavx2, wexx,ua2_(niw),freqw1,q_r(3),qk(3)
      logical,parameter :: debug=.false.,timemix=.false.
      logical ::   oncew, onceww, eibz4sig  
      real(8),allocatable:: we_(:,:),wfac_(:,:)
      complex(8),allocatable:: w3p(:),wtff(:)
      logical:: tote=.false., hermitianW
      real(8),allocatable:: vcoud_(:),wfft(:)
      logical:: iprx,cmdopt0
      integer:: ixx,ixc,icount
      real(8),parameter:: pi=4d0*datan(1d0), fpi=4d0*pi, tpi=8d0*datan(1d0),ddw=10d0
      integer:: ncount,kxold,nccc
      integer,allocatable:: ispc(:),kxc(:),irotc(:),ipc(:),krc(:)
      integer,allocatable:: nwxic(:), nwxc(:), nt_maxc(:),irkip(:,:,:,:)
c      if(verbose()>=90) debug=.true.
!!----------------------------------------------------------------
!!!!!!!!!!!
      if(npm==2) call rx('sxcf: npm=2 need to be examined')
!!!!!!!!!!!

      allocate(zsecall(ntq,ntq,nq,nspinmx)) !, coh(ntq,nq) ) kount(nqibz,nq),
      zsecall = 0d0
      if(ixc==3.and.nctot==0) return
      
!  We divide irkip_all into irkip for nodes. irkip is dependent on rank.
!  Total number of none zero irkip for all ranks is the number of nonzero irkip_all
      allocate(irkip(nspinmx,nqibz,ngrp,nq)) ! nrkip is weight correspoinding to irkip for a node.
      call MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nq,  irkip)
!  icount mechanism
      IcountBlock: Block
      ncount=count(irkip/=0)
      allocate(ispc(ncount),kxc(ncount),irotc(ncount),ipc(ncount),krc(ncount))
      allocate(nwxic(ncount), nwxc(ncount), nt_maxc(ncount))
      iqini = 1
      iqend = nqibz             !no sum for offset-Gamma points.
      icount=0
      do 1300 kx = iqini,iqend 
      do 1200 isp = 1,nspinmx !empty run to get index for icount ordering
         if(sum(irkip(isp,kx,:,:))==0) cycle ! next kx
      do 1400 irot = 1,ngrp !over rotations irot ===
         if(sum(irkip(isp,kx,irot,:))==0) cycle ! next ip
      do 1500 ip = 1,nq         
         kr = irkip(isp,kx,irot,ip) ! index for rotated kr in the FBZ
         if(kr==0) cycle
         icount=icount+1  
         ispc(icount)=isp !icount specify isp,kx,irot,iq. (kx,irot) gives kr in the all FZ.
         kxc(icount)=kx
         irotc(icount)=irot
         ipc(icount)=ip
         krc(icount)=kr
         qibz_k = qibz(:,kx)
         q(1:3)= qip(1:3,ip)
         eq = readeval(q,isp)
         omega(:) = eq(itq(:))  !1:ntq
         qbz_kr= qbz (:,kr)     !rotated qbz vector. 
         qk =  q - qbz_kr        
         ekq = readeval(qk, isp) 
         ekc(nctot+1:nctot+nband) = ekq (1:nband)
         nt0 = count(ekc<ef) 
         nt0p = count(ekq<ef+ddw*esmr) +nctot 
         nt0m = count(ekq<ef-ddw*esmr) +nctot
         ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
         write(6,*) icount, ispc(icount),kxc(icount),' irot ',irot,ip,kr
         if(exchange) cycle
!   correlation case only nwxi,nwx nt_max
         ebmx=1d10              !this is needed probably because to fill 1d99 for ekc(i) above boundary.
         nbmxe = count(ekc<ebmx)-nctot !nocc (ekc,ebmx,nstatetot)-nctot!
         nbmax  = min(nband,nbmxe) 
         nstate = nctot + nbmax ! = nstate for the case of correlation
         call get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,
     i        nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,
     o        nwxi,nwx,nt_max)  ! Get index nwxi nwx nt_max. 
         nwxic(icount)=nwxi     ! get_nwx is not written clearly, but works and not time-consuming.
         nwxc(icount)=nwx
         nt_maxc(icount)=nt_max
 1500 continue
 1400 continue
 1200 continue
 1300 continue
      if(icount/=count(irkip/=0)) call rx('sxcf: icount/=count(irkip/=0)')
      EndBlock IcountBlock
      
! open WV* files     
      if(.not.exchange) then
         allocate(ifrcw(iqini:iqend),ifrcwi(iqini:iqend))
         do kx=iqini,iqend
           if(any(kx==kxc)) then !only for requied files for this rank
           open(newunit=ifrcw(kx),  file='WVR.'//charnum5(kx), action='read',form='unformatted',
     &           status='unknown',access='direct',recl=mrecl)
           open(newunit=ifrcwi(kx), file='WVI.'//charnum5(kx), action='read',form='unformatted',
     &           status='unknown',access='direct',recl=mrecl)
           endif
         enddo
      endif
      
!  main loop do 3030    
      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot,isp) ! core
      kxold=-9999
      nccc= max(ncount/50,1) !just for printing
      do 3030 icount=1,ncount   !we only consider bzcase()==1 now
         isp =ispc(icount)
         kx  =kxc(icount) !for W(kx), kx is irreducible
         irot=irotc(icount)
         ip  =ipc(icount)
         kr  =krc(icount) !=irkip(isp,kx,irot,ip) runs all the k mesh points required for q(ip).
         nt_max=nt_maxc(icount)
         nwxi=nwxic(icount)
         nwx =nwxc(icount)
         if(mod(icount,nccc)==0) 
     &     write(6,"('sxcf: isp kx ip',3i5,' icount of ncount=',2i5)") isp,kx,ip,icount,ncount
         zsec => zsecall(:,:,ip,isp)
         qibz_k = qibz(:,kx)
         q(1:3)= qip(1:3,ip)
         eq = readeval(q,isp)
         omega(:) = eq(itq(:))  !1:ntq
         qbz_kr= qbz (:,kr)     !rotated qbz vector. 
         qk =  q - qbz_kr        
         ekq = readeval(qk, isp) 
         ekc(nctot+1:nctot+nband) = ekq (1:nband)
         nt0 = count(ekc<ef) 
         nt0p = count(ekq<ef+ddw*esmr) +nctot 
         nt0m = count(ekq<ef-ddw*esmr) +nctot
         ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
         wtt = wk(kr)
         if(eibz4sig()) wtt=wtt*nrkip(isp,kx,irot,ip)
c         if(debug) write(6,*)' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
         if(kxold/=kx) then
            call Readvcoud(qibz_k,kx,NoVcou=.false.) !Readin ngc,ngb,vcousq ! Coulomb matrix
            call Setppovlz(qibz_k,matz=.true.) !ppovlz overlap matrix
            if(debug) write(6,*) ' sxcf_fal2sc: ngb ngc nbloch=',ngb,ngc,nbloch
            kxold =kx
         endif
!! Get zmel(ib,itpp,it) = <M(qbz_kr,ib) phi(itpp,q-qbz_kr) |phi(q(ip),it)> , qbz_kr= irot(qibz_k)
         if(exchange) then
            nbmax = nt0p-nctot
         else
            ebmx=1d10  !this is needed probably because to fill 1d99 for ekc(i) above boundary.
            nbmxe = count(ekc<ebmx)-nctot !nocc (ekc,ebmx,nstatetot)-nctot!
            nbmax  = min(nband,nbmxe) 
         endif
         nstate = nctot + nbmax ! = nstate for the case of correlation
         call Get_zmel_init(q,qibz_k,irot,qbz_kr,isp, nbmax,ntqxx,nctot,ncc=0,iprx=debug)
!---------
        ExchangeSelfEnergy: Block
         real(8):: wfacx
         if(exchange) then      
            allocate(vcoud_(ngb),wtff(nctot+nbmax),w3p(nctot+nbmax)) 
            vcoud_= vcoud
            if(kx == iqini) vcoud_(1) = wklm(1)* fpi*sqrt(fpi) /wk(kx)
                    !voud_(1) is effective v(q=0) averaged in the Gamma cell.
            wtff= [(wfacx(-1d99, ef, ekc(it), esmr),it=1,nctot+nbmax)]
            do itpp= 1,ntqxx
            do itp = 1,ntqxx !S[j=1,nbloch]  z1p(j,t,n) <B(rk,j) psi(q-rk,n) |psi(q,t')>
               w3p = [(sum(dconjg(zmel(:,it,itp)) * vcoud_(:)* zmel(:,it,itpp)),it=1,nctot+nbmax)]
               w3p(nctot+1:nctot+nbmax) = w3p(nctot+1:nctot+nbmax) * wtff(nctot+1:nctot+nbmax)
               if(corehole) w3p(1:nctot) = w3p(1:nctot) * wcorehole(1:nctot,isp) 
               zsec(itp,itpp) = zsec(itp,itpp) - wtt * sum( w3p(:) )
            enddo
            enddo
            deallocate(vcoud_,wtff,w3p)
            cycle               
         endif                  ! end of if(exchange)
        EndBlock ExchangeSelfEnergy
         
         if(timemix) call timeshow("33333 k-cycle")
         
!     ! Integration along imag axis for zwz(omega) for given it,itp,itpp
!     ! itp  : left-hand end of expternal band index.
!     ! itpp : right-hand end of expternal band index.
!     ! it   : intermediate state of G.
!     !===  See Eq.(55) around of PRB76,165106 (2007)
!
!---------
        CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
         real(8):: esmrx(nstate),omegat(ntqxx),wgtim(0:npm*niw,ntqxx,nstate)
         esmrx(1:nctot)=0d0
         esmrx(nctot+1:nstate)=esmr
         omegat(1:ntqxx) = omega(1:ntqxx)
         do itp = 1,ntqxx
         do  it = 1,nstate
            we = .5d0*(omegat(itp)-ekc(it))
            call wintzsg_npm_wgtim(npm, ua_,expa_, we,esmrx(it), wgtim(:,itp,it))
         enddo   !Integration weight wgtim along im axis for zwz(0:niw*npm) 
         enddo
         allocate( zwz(nstate,ntqxx,ntqxx), zw (nblochpmx,nblochpmx))
         do ixx=0,niw !niw is ~10. ixx=0 is for omega=0 nw_i=0 (Time reversal) or nw_i =-nw
            if(ixx==0) then ! at omega=0 ! nw_i=0 (Time reversal) or nw_i =-nw
               read(ifrcw(kx),rec=1+(0-nw_i)) zw ! direct access read Wc(0) = W(0) - v
            elseif(ixx>0) then ! 
               read(ifrcwi(kx),rec=ixx) zw ! direct access read Wc(i*omega)=W(i*omega)-v
            endif
            call matzwz2(2,zw(1:ngb,1:ngb),zmel,ntqxx,nstate,ngb,zwz) !zwz = zmel*(W(i w)-v)*zmel
            do istate = 1,nstate
               Block
               complex(8):: zwx(ntqxx,ntqxx)
               zwx = zwz(istate,1:ntqxx,1:ntqxx)
               zwz(istate,:,:)=(zwx+dconjg(transpose(zwx)))/2d0 ! w(iw) + w(-iw) Harmitian part
c     if(npm==2)zwz(ix+niw,istate,:,:)=(zwx-dconjg(transpose(zwx)))/2d0/img!w(iw)-w(-iw) Anti
               EndBlock
            enddo
            do itp = 1,ntqxx
              zsec(itp,:)=zsec(itp,:)+ wtt*matmul(wgtim(ixx,itp,:),zwz(:,itp,:)) !itp,itpp
            enddo
         enddo
         deallocate(zw)
         deallocate(zwz)
        EndBlock CorrelationSelfEnergyImagAxis
!---------
        CorrelationSelfEnergyRealAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
         real(8):: we_(nt_max,ntqxx),wfac_(nt_max,ntqxx)
         complex(8):: zw3(ngb,ngb,0:2)
         integer:: ixss(nt_max,ntqxx),iirx(ntqxx)
         logical:: ititpskip(nt_max,ntqxx)
         complex(8):: zw(nblochpmx,nblochpmx),img=(0d0,1d0)
         if(debug)write(6,*)' CorrelationSelfEnergyRealAxis: Block '
         call weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,
     i        ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt,
     o        we_,wfac_,ixss,ititpskip,iirx)
         Block
         real(8):: wgt3(0:2,nt_max,ntqxx)    ! 3-point interpolation weight for we_(it,itp) 
         complex(8)::zadd(ntqxx),zw33(ngb,ngb) ! ixss is starting index of omega
         integer:: iwgt3(nt_max,ntqxx),i1,i2,iw,ikeep,ix
         integer:: ncountp(nwxi:nwx),icountp,ncoumx
         integer,allocatable:: itc(:,:),itpc(:,:)
         real(8),allocatable:: wgt3p(:,:,:)
         iwgt3=0
         do itp=1,ntqxx
         do it=1,nt_max
            if(ititpskip(it,itp)) cycle
            ixs= ixss(it,itp) !we_ is \omega_\epsilon in Eq.(55).
            if(npm==1) call alagr3z2wgt(we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
            if(npm==2) call alagr3zwgt (we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
            iwgt3(it,itp)=iirx(itp)*(ixs+1-2) !starting omega index ix for it,itp
         enddo                                !iirx=1 for npm=1 I think
            if(iirx(itp)/=1) call rx('sxcf: iirx=-1(TR breaking) is not yet implemented')
         enddo
!! icount mechanism for sum
!  We take sum for ix,it,itp whereas W(we_(it,itp))=\sum_{i=0}^2 W(:,:,ix+i)*wgt3(i)         
         do ix = nwxi,nwx 
           icountp=0
           do itp=1,ntqxx
           do it=1,nt_max
            if(ititpskip(it,itp)) cycle
            if(iwgt3(it,itp)/=ix) cycle
            icountp=icountp+1
           enddo
           enddo
           ncountp(ix)=icountp
         enddo
         ncoumx=maxval(ncountp)
         allocate(itc(ncoumx,nwxi:nwx),itpc(ncoumx,nwxi:nwx),wgt3p(0:2,ncoumx,nwxi:nwx))
         do ix = nwxi,nwx 
           icountp=0
           do itp=1,ntqxx
           do it=1,nt_max
            if(ititpskip(it,itp)) cycle
            if(iwgt3(it,itp)/=ix) cycle
            icountp=icountp+1
            itc (icountp,ix)=it
            itpc(icountp,ix)=itp
            wgt3p(:,icountp,ix)=wgt3(:,it,itp)
           enddo
           enddo
         enddo   
!   ix-shift mechanism for reading zw(:,:,ix)
         ikeep=99999
         do ix = nwxi,nwx       !Set zw3(:,:,0:2) is for ix,ix+1,ix+2
           if(.not.(any(iwgt3(:,:)==ix))) cycle
           do iw=0,2
              if(ikeep+1==ix.and.iw<2) then
                 zw3(:,:,iw)=zw3(:,:,iw+1) 
              elseif(ikeep+2==ix.and.iw<1) then
                 zw3(:,:,iw)=zw3(:,:,iw+2)
              else
                 read(ifrcw(kx),rec=iw+ix-nw_i+1) zw ! direct access Wc(omega) = W(omega) - v
                 if(hermitianW) zw3(:,:,iw)=(zw(1:ngb,1:ngb)+transpose(dconjg(zw(1:ngb,1:ngb))))/2d0
                 if(.not.hermitianW) zw3(:,:,iw)=zw(1:ngb,1:ngb)
              endif   
           enddo
           ikeep=ix
!         ! zw3 should contain zw(ix+0:ix+2)
           do icountp=1,ncountp(ix)
              itp=itpc(icountp,ix)
              it =itc(icountp,ix)  !zw33 gives interpolated value of W(we_(it,itp))
              zw33 = zw3(:,:,0)*wgt3p(0,icountp,ix) 
     &             + zw3(:,:,1)*wgt3p(1,icountp,ix) 
     &             + zw3(:,:,2)*wgt3p(2,icountp,ix) 
              zadd(:) = matmul(matmul(dconjg(zmel(:,it,itp)),zw33),zmel(:,it,:)) !itpp
              zadd(itp)= dreal(zadd(itp))+img*min(dimag(zadd(itp)),0d0) !forced to be negative im part
              zsec(itp,:)=zsec(itp,:)+ zadd(:)* wfac_(it,itp) !itpp
            enddo
         enddo
         Endblock         
        EndBlock CorrelationSelfEnergyRealAxis
!         
         call Deallocate_zmel()
         if(verbose()>50) call timeshow("11after alagr3z iw,itp,it cycles")
 3030 continue
      
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
      integer:: itini,iii,it,itend,wp,ixs=-9999,itp,iwp,nt0p
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
c              write(6,*)'ixs+1 nwx=',ixs+1,nwx
              call rx( ' sxcf: ixs+1>nwx yyy2a')
            endif
            if(omg < ef .and. abs(ixs+1)> abs(nwxi) ) then
c              write(6,*)'ixs+1 nwxi=',ixs+1,nwxi
              call rx( ' sxcf: ixs-1<nwi yyy2b')
            endif
          endif
        enddo     
      enddo
      end subroutine weightset4intreal
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      subroutine get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,
     i nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,
     o nwxi,nwx,nt_max)
      implicit none
      intent(in)::       omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,
     i     nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug
      intent(out)::
     o nwxi,nwx,nt_max
!> Determine indexes of a range for calculation. !! It is better to clean this up...
      integer:: nctot,nw_i,nw,nstate,nt0p,nt0m,ntq, nband,ntqxx
      real(8):: esmr,ef,ekc(nctot+nband),wfaccut,freq_r(nw_i:nw)
      real(8):: wfac,wfacx2,we,weavx2,esmrx,wexx
      real(8),pointer::omg
      real(8),target:: omega(ntq)
      integer:: nt_max,nwxi,nwx,itp,it,itini,itend,iwp,ixs=-9999,ixsmin,ixsmx,verbose
      logical::debug
!!     maximum ixs reqired.
      ixsmx =0
      ixsmin=0
      do 301 itp = 1,ntqxx
        omg => omega(itp) 
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
            write (*,*)'nw_i ixsmin wexx',nw_i,ixsmin,wexx,' omg ekc(it) ef ', omg,ekc(it),ef
            call rx( ' sxcf 222: |w-e| out of range')
          endif
 311    continue
 301  continue               
      if(nw_i==0) then          !time reversal
        nwxi = 0
        nwx  = max(ixsmx+1,ixsmin+1)
      else                      !no time revarsal working?
        nwxi = -ixsmin-1
        nwx  =  ixsmx+1
      endif
      if(nwx > nw .or. nwxi < nw_i ) call rx( ' get_nwx : |w-e| > max(w)')
c      if(verbose()>50)call timeshow("10before alagr3z iw,itp,it ")
      nt_max=nt0p               !initial nt_max
      do 401 itp = 1,ntqxx
        omg => omega(itp)
        if (omg > ef) then
          do  it = nt0m+1,nstate ! nt0m corresponds to efm
            wfac = wfacx2 (ef,omg, ekc(it),esmr)
            if(wfac>wfaccut) then
              if (it > nt_max) nt_max=it ! nt_max is  unocc. state
            endif               ! that ekc(it>nt_max)-omega > 0
          enddo                 ! so it > nt_max does not contribute to omega pole integral
        endif
 401  continue               
      end subroutine get_nwx

      end module m_sxcfsc
