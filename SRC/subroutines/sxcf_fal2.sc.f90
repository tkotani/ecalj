module m_sxcfsc !self-energy calculation 
  use m_readqg,only   : Readqg0
  use m_readeigen,only: Readeval
  use m_keyvalue,only   : Getkeyvalue
  use m_zmel,only : Get_zmel_init, Setppovlz, Setppovlz_chipm, Deallocate_zmel, get_zmel_modex0, zmel
  use m_itq,only: itq,ntq
  use m_genallcf_v3,only: nlmto,nspin,nctot,niw,ecore !,symgg
  use m_read_bzdata,only: qibz,qbz,wk=>wbz,nqibz,nqbz,wklm,lxklm,nq0i, wqt=>wt,q0i, irk
  use m_readVcoud,only:   Readvcoud, vcoud,ngb,ngc
  use m_readfreq_r,only: freq_r, nw_i,nw,freqx,wx=>wwx,nblochpmx,mrecl,expa_,npm,nprecx
  use m_rdpp,only: Rdpp, nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
  use m_readqg,only:  ngpmx,ngcmx
  use m_readhbe,only: nband,mrecg
  use m_hamindex,only: ngrp
  !  use m_eibzhs,only: nrkip=>nrkip_all,irkip_all
!  use m_read_bzdata,only: qbz,nqibz,ginv,irk,nqbz
  use m_readgwinput,only: ua_,  corehole,wcorehole
  use m_mpi,only: MPI__sxcf_rankdivider
  use m_ftox
  implicit none
  !---------------------------------      
  public sxcf_scz, zsecall
  complex(8),allocatable,target:: zsecall(:,:,:,:) !output
  private
contains

  subroutine sxcf_scz(qvec,ef,esmr,nq,exchange,nbandmx,ixc,nspinmx) !,jobsw
    intent(in)        qvec,ef,esmr,nq,exchange,nbandmx,ixc,nspinmx !,jobsw
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
    !!xxx jobsw switch. We now support only mode=3-----------------
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
    !!     q     =qvec(:,iq)  = q-vector in SEc(q,t). 
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
    real(8) :: ef,esmr, qvec(3,nq)
    real(8):: ebmx
    complex(8),pointer::zsec(:,:)
    complex(8),pointer::ww(:,:)
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
    real(8) :: wtt,wfac,we!,esmrx
    real(8) :: qvv(3),eq(nband),omega(ntq),quu(3),freqw,ratio
    real(8) :: qibz_k(3),qbz_kr(3),vc,omega0,omg
    complex(8),allocatable,target:: zwz(:,:,:),zw(:,:)
    real(8), parameter :: wfaccut=1d-8,tolq=1d-8
    complex(8), parameter :: img=(0d0,1d0)
    character(10) :: i2char
    real(8)::polinta, wfacx, wfacx2, weavx2, wexx,ua2_(niw),freqw1,q_r(3),qk(3)
    logical,parameter :: debug=.false.,timemix=.true.
    logical ::   oncew, onceww !, eibz4sig  
    real(8),allocatable:: we_(:,:),wfac_(:,:)
    complex(8),allocatable:: w3p(:),wtff(:)
    logical:: tote=.false.!, hermitianW
    real(8),allocatable:: vcoud_(:),wfft(:)
    logical:: iprx,cmdopt0
    integer:: ixx,ixc,icount,ndivmx,ns1,ns2,ns1c,ns2c,ns1v
    real(8),parameter:: pi=4d0*datan(1d0), fpi=4d0*pi, tpi=8d0*datan(1d0),ddw=10d0
    integer:: ncount,kxold,nccc,icount0
    integer,allocatable:: ispc(:),kxc(:),irotc(:),ipc(:),krc(:),nstateMax(:),nstti(:),nstte(:)
    integer,allocatable:: nwxic(:), nwxc(:), nt_maxc(:),irkip(:,:,:,:)
    complex(8),allocatable:: zmelc(:,:,:)
    integer,allocatable::ndiv(:),nstatei(:,:),nstatee(:,:)
    !!----------------------------------------------------------------
    if(npm==2) call rx('sxcf: npm=2 need to be examined')
    allocate(zsecall(ntq,ntq,nq,nspinmx)) !, coh(ntq,nq) ) kount(nqibz,nq),
    zsecall = 0d0
    if(ixc==3.and.nctot==0) return
    !  We divide irkip_all into irkip for nodes. irkip is dependent on rank.
    !  Total number of none zero irkip for all ranks is the number of nonzero irkip_all
    block
      integer:: irkip_all(nspinmx,nqibz,ngrp,nq),iqq,is
      do is = 1,nspinmx
         do iqq=1,nq
            irkip_all(is,:,:,iqq)=irk
         enddo
      enddo
      allocate(    irkip(nspinmx,nqibz,ngrp,nq)) ! nrkip is weight correspoinding to irkip for a node.
      call MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nq,  irkip)
    endblock
!
    PreIcountBlock: Block!Get nstateMax(ncount),ndiv(icount),nstatei(j,icount),nstatee(j,icount)
      integer:: ndivide,nstateavl,nnn,nloadav,nrem,idiv,j
      integer,allocatable:: nload(:)
      ncount=count(irkip/=0)
      allocate(nstateMax(ncount))
      iqini = 1
      iqend = nqibz             
      icount=0
      kxloop:do kx = iqini,iqend !quick loop
         isploop: do isp = 1,nspinmx !empty run to get index for icount ordering
            if(sum(irkip(isp,kx,:,:))==0) cycle ! next kx
            irotloop: do irot = 1,ngrp !over rotations irot ===
               if(sum(irkip(isp,kx,irot,:))==0) cycle ! next ip
               iqloop: do 1150 ip = 1,nq         
                  kr = irkip(isp,kx,irot,ip) ! index for rotated kr in the FBZ
                  if(kr==0) cycle
                  icount=icount+1
                  q(1:3)= qvec(1:3,ip)
                  qbz_kr= qbz (:,kr)     !rotated qbz vector. 
                  qk =  q - qbz_kr        
                  ekq = readeval(qk, isp) 
                  ekc(nctot+1:nctot+nband) = ekq (1:nband)
                  nt0p = count(ekq<ef+ddw*esmr) +nctot 
                  if(exchange) then
                     nstateMax(icount) = nt0p
                  else   
                     ebmx=1d10 !this is needed probably for filling 1d99 for ekc(i) above boundary.
                     nbmxe = count(ekc<ebmx)-nctot !nocc (ekc,ebmx,nstatetot)-nctot!
                     nbmax  = min(nband,nbmxe) 
                     nstateMax(icount) = nctot + nbmax ! = nstate for the case of correlation
                  endif
1150           enddo iqloop
            enddo irotloop
         enddo isploop
      enddo kxloop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
                      ! ndivide = 2 
      nstateavl = 16  ! nstateavl=max(sum(nstatemax)/(ncount*ndivide),1)
      if(ixc==3) nstateavl= maxval(nstatemax)
      ! size of average load of middle states (in G)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(ndiv(ncount))
      ndiv = (nstatemax-1)/nstateavl + 1  !number of division for icount
      ndivmx = maxval(ndiv)
      allocate(nstatei(ndivmx,ncount),nstatee(ndivmx,ncount),nload(ndivmx))
      do icount=1,ncount
         nnn = nstatemax(icount)    !total number of middle states for given icount
         ndiv(icount) = (nnn-1)/nstateavl + 1  !number of division for icount
         nloadav = nnn/ndiv(icount) !number of average load of middle states
         nrem = nnn - ndiv(icount)*nloadav !remnant count
         nload(1:nrem) = nloadav+1
         nload(nrem+1:ndiv(icount)) = nloadav !write(6,ftox)'nload=',nload(1:ndiv(icount))
         nstatei(:,icount)= [(sum(nload(1:idiv-1))+1,idiv=1,ndiv(icount))]!init index for(idiv,icount)
         nstatee(:,icount)= [(sum(nload(1:idiv)),    idiv=1,ndiv(icount))]!end  index
      enddo
!      do icount=1,ncount
!      do j=1,ndiv(icount)
!        write(6,ftox)'nnnx icou ndiv=',icount,j,nstatei(j,icount),nstatee(j,icount),nstatemax(icount)
!      enddo
!      enddo
      deallocate(nload,nstatemax)
    EndBlock PreIcountBlock
    write(6,*)'nnn init ncount=',ncount
    
    ncount = ncount*ndivmx
    ! icount mechanism for sum in MAINicountloop 3030
    IcountBlock: Block !quick loop to gather index sets for main loop
      integer:: idiv
!      ncount=count(irkip/=0)
      allocate(ispc(ncount),kxc(ncount),irotc(ncount),ipc(ncount),krc(ncount))
      allocate(nwxic(ncount), nwxc(ncount), nt_maxc(ncount),nstateMax(ncount))
      allocate(nstti(ncount),nstte(ncount))
      iqini = 1
      iqend = nqibz             !no sum for offset-Gamma points.
      icount=0
      icount0=0
      do 130 kx = iqini,iqend !this is empty run to get index for icount ordering
         do 120 isp = 1,nspinmx 
            if(sum(irkip(isp,kx,:,:))==0) cycle ! next kx
            do 140 irot = 1,ngrp !over rotations irot ===
               if(sum(irkip(isp,kx,irot,:))==0) cycle ! next ip
               do 150 ip = 1,nq         
                  kr = irkip(isp,kx,irot,ip) ! index for rotated kr in the FBZ
                  if(kr==0) cycle
                  icount0=icount0+1
                  do idiv=1,ndiv(icount0) !icount loop have further division by ndiv
                     icount=icount+1
                     nstti(icount)=nstatei(idiv,icount0) ![nstti,nstte] specify range of int. states.
                     nstte(icount)=nstatee(idiv,icount0) !
                     ispc(icount)=isp!icount specify isp,kx,irot,iq. (kx,irot) gives kr in the all FZ.
                     kxc(icount)=kx
                     irotc(icount)=irot
                     ipc(icount)=ip
                     krc(icount)=kr
                     qibz_k = qibz(:,kx)
                     q(1:3)= qvec(1:3,ip)
                     eq = readeval(q,isp)
                     omega(:) = eq(itq(:))  !1:ntq
                     qbz_kr= qbz (:,kr)     !rotated qbz vector. 
                     qk =  q - qbz_kr        
                     ekq = readeval(qk, isp) 
                     ekc(nctot+1:nctot+nband) = ekq (1:nband)
                     nt0  = count(ekc<ef) 
                     nt0p = count(ekq<ef+ddw*esmr) +nctot 
                     nt0m = count(ekq<ef-ddw*esmr) +nctot
                     ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
                     !write(6,*) icount, ispc(icount),kxc(icount),' irot ',irot,ip,kr
                     if(exchange) then
                        nstateMax(icount) = nt0p
                     else   
                        ebmx=1d10 !this is because filling 1d99 for ekc(i) above boundary.
                        nbmxe = count(ekc<ebmx)-nctot !nocc (ekc,ebmx,nstatetot)-nctot!
                        nbmax  = min(nband,nbmxe) 
                        nstateMax(icount) = nctot + nbmax ! = nstate for the case of correlation
                        call get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstateMax(icount),freq_r,&
                             nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,nwxi,nwx,nt_max)
                        !Get index nwxi nwx nt_max.
                        ! get_nwx is not written clearly, but works and not time-consuming.
                        nwxic(icount)=nwxi 
                        nwxc(icount)=nwx
                        nt_maxc(icount)=nt_max
                     endif
                  enddo
150            enddo
140         enddo
120      enddo
130   enddo
      if(icount0/=count(irkip/=0)) call rx('sxcf: icount/=count(irkip/=0)')
      ncount=icount
    EndBlock IcountBlock
    write(6,*)'nnn dev  ncount=',ncount
    
    if(.not.exchange) then! Read WV* containing W-v in MPB
       allocate(ifrcw(iqini:iqend),ifrcwi(iqini:iqend))
       do kx=iqini,iqend
          if(any(kx==kxc)) then !only for requied files for this rank
             open(newunit=ifrcw(kx), file='WVR.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
             open(newunit=ifrcwi(kx),file='WVI.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
          endif
       enddo
    endif
    
    if(nctot/=0) ekc(1:nctot)= ecore(1:nctot,isp) ! core
    kxold=-9999
    !nccc= max(ncount/50,1) !just for printing
    MAINicountloop: do 3030 icount=1,ncount   !we only consider bzcase()==1 now
       !write(6,*)'do 3030 icount=',icount
       isp =ispc(icount)
       kx  =kxc(icount) !for W(kx), kx is irreducible
       irot=irotc(icount)
       ip  =ipc(icount)
       kr  =krc(icount) !=irkip(isp,kx,irot,ip) runs all the k mesh points required for q(ip).
       nt_max=nt_maxc(icount)
       nwxi=nwxic(icount)
       nwx =nwxc(icount)
       qibz_k = qibz(:,kx)
       q(1:3)= qvec(1:3,ip)
       eq = readeval(q,isp)
       omega(:) = eq(itq(:))  !1:ntq
       qbz_kr= qbz (:,kr)     !rotated qbz vector. 
       qk =  q - qbz_kr        
       ekq = readeval(qk, isp) 
       ekc(nctot+1:nctot+nband) = ekq (1:nband)
       nt0  = count(ekc<ef) 
       nt0p = count(ekq<ef+ddw*esmr) +nctot 
       nt0m = count(ekq<ef-ddw*esmr) +nctot
       ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
       zsec => zsecall(1:ntqxx,1:ntqxx,ip,isp)
       wtt = wk(kr)
       !if(eibz4sig()) wtt=wtt*nrkip(isp,kx,irot,ip)
       !! Get zmel(ib,itpp,it) = <M(qbz_kr,ib) phi(itpp,q-qbz_kr) |phi(q(ip),it)> , qbz_kr= irot(qibz_k)
       if(kxold/=kx) then
          call Readvcoud(qibz_k,kx,NoVcou=.false.) !Readin ngc,ngb,vcoud ! Coulomb matrix
          call Setppovlz(qibz_k,matz=.true.) !Set ppovlz overlap matrix in m_zmel
          if(debug) write(6,*) ' sxcf_fal2sc: ngb ngc nbloch=',ngb,ngc,nbloch
          kxold =kx
       endif
       nstate = nstateMAX(icount) ! for the case of correlation
       ns1=nstti(icount) ! Range of middle states is [ns1:ns2] 
       ns2=nstte(icount) ! nstte(icount) is upper bound of middle states
       write(6,ftox)'do3030:isp kx ip',isp,kx,ip,'icou/ncou=',icount,ncount,'ns1:ns2=',ns1,ns2
       call get_zmel_modex0(ns1,1,isp,isp) !set lower bound of middle state
       call Get_zmel_init(q,qibz_k,irot,qbz_kr,isp, ns2-nctot,ntqxx,nctot,ncc=0,iprx=debug)!Get zmel
       Exchangemode: if(exchange) then      
          ExchangeSelfEnergy: Block
            real(8):: wfacx
            allocate(vcoud_(ngb),wtff(ns1:ns2),w3p(ns1:ns2)) !range of middle states ns1:ns2
            vcoud_= vcoud
            if(kx == iqini) vcoud_(1) = wklm(1)* fpi*sqrt(fpi) /wk(kx) !voud_(1) is effective v(q=0) averag in the Gamma cell.
            wtff= [(wfacx(-1d99, ef, ekc(it), esmr),it=ns1,ns2)]
            do itpp= lbound(zsec,2),ubound(zsec,2)
               do itp = lbound(zsec,1),ubound(zsec,1)
                  w3p(ns1:ns2)=[(sum(dconjg(zmel(:,it,itp))*vcoud_(:)*zmel(:,it,itpp)),it=1,ns2-ns1+1)]
                  ns1v= max(nctot+1,ns1) ! minimum index of valence for core+valence
                  w3p(ns1v:ns2) = w3p(ns1v:ns2) * wtff(ns1v:ns2)
                  if(corehole) then !not checked well
                     ns2c= min(ns2,nctot) !ns2c upper limit of core index of core+valence
                     ns1c= ns1 !ns1c lower limit of core index of core+valence
                     w3p(ns1c:ns2c) = w3p(ns1c:ns2c) * wcorehole(ns1c:ns2c,isp)
                  endif
                  zsec(itp,itpp) = zsec(itp,itpp) - wtt * sum( w3p(:) )
               enddo
            enddo
            deallocate(vcoud_,wtff,w3p)
            if(timemix) call timeshow("ExchangeSelfEnergy cycle")
          EndBlock ExchangeSelfEnergy
          cycle  !end of exchange mode             
       endif Exchangemode
       !     ! Integration along imag axis for zwz(omega) for given it,itp,itpp
       !     ! itp  : left-hand end of expternal band index.
       !     ! itpp : right-hand end of expternal band index.
       !     ! it   : intermediate state of G.
       !     !===  See Eq.(55) around of PRB76,165106 (2007)
       !
       zmelcww: Block !range of middle states is [ns1:ns2] instead of [1:nstate]. 2022-8-28
         complex(8):: zmelcww(1:ntqxx,ns1:ns2,1:ngb)
         zmelcww=0d0
         allocate(zmelc(1:ntqxx,ns1:ns2,1:ngb)) !1:nstate,1:ngb))
         forall(itp=1:ntqxx) zmelc(itp,:,:)=transpose(dconjg(zmel(:,:,itp))) 
         CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
           real(8):: esmrx(nstate), omegat(ntqxx), wgtim(0:npm*niw,ntqxx,ns1:ns2)
           complex(8),target:: zw(nblochpmx,nblochpmx),zwz(ns1:ns2,ntqxx,ntqxx)
           logical:: init=.true.
           if(timemix) call timeshow(" CorrelationSelfEnergyImagAxis:")
           ns2c= min(ns2,nctot)
           ns1c= max(min(ns1,nctot),1)
           ns1v= max(nctot+1,ns1) !    !print *,'ns1c,ns2c,ns1v=',ns1c,ns2c,ns1v
           esmrx(ns1c:ns2c)= 0d0
           esmrx(ns1v:ns2) = esmr
           omegat(1:ntqxx) = omega(1:ntqxx)
           itpdo:do itp = lbound(zsec,1), ubound(zsec,1)
              do     it = lbound(zmelc,2),ubound(zmelc,2)
                 we = .5d0*(omegat(itp)-ekc(it))
                 call wintzsg_npm_wgtim(npm, ua_,expa_, we,esmrx(it), wgtim(:,itp,it))
              enddo   !Integration weight wgtim along im axis for zwz(0:niw*npm) 
           enddo itpdo
           iwimag:do ixx=0,niw !niw is ~10. ixx=0 is for omega=0 nw_i=0 (Time reversal) or nw_i =-nw
              if(ixx==0) then ! at omega=0 ! nw_i=0 (Time reversal) or nw_i =-nw
                 read(ifrcw(kx),rec=1+(0-nw_i)) zw ! direct access read Wc(0) = W(0) - v
              elseif(ixx>0) then ! 
                 read(ifrcwi(kx),rec=ixx) zw ! direct access read Wc(i*omega)=W(i*omega)-v
              endif
              ww=>zw(1:ngb,1:ngb)
              call matmaw(zmelc,ww,zmelcww, size(zmelc,1)*size(zmelc,2), size(ww,1), size(ww,2),&
                   &           wtt*wgtim(ixx,:,:))
           enddo iwimag
         EndBlock CorrelationSelfEnergyImagAxis 
         CorrelationSelfEnergyRealAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
           real(8):: we_(nt_max,ntqxx),wfac_(nt_max,ntqxx)
           integer:: ixss(nt_max,ntqxx),iirx(ntqxx)
           logical:: ititpskip(nt_max,ntqxx)
           complex(8),target:: zw(nblochpmx,nblochpmx)
           if(timemix) call timeshow(" CorrelationSelfEnergyRealAxis:")
           if(debug)write(6,*)' CorrelationSelfEnergyRealAxis: Block '
           call weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,&
                ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt,&
                we_,wfac_,ixss,ititpskip,iirx)
           CorrR2:Block
             real(8):: wgt3(0:2,nt_max,ntqxx)    ! 3-point interpolation weight for we_(it,itp) 
             complex(8)::zadd(ntqxx),wv33(ngb,ngb) ! ixss is starting index of omega
             complex(8):: wv3(ngb,ngb,0:2)
             integer:: iwgt3(nt_max,ntqxx),i1,i2,iw,ikeep,ix
             integer:: nit_(ntqxx,nwxi:nwx),icountp,ncoumx,iit
             integer,allocatable:: itc(:,:,:),itpc(:,:)
             real(8),allocatable:: wgt3p(:,:,:)
             if(timemix) call timeshow(" CorrR2:")
             iwgt3=0
             do itp=1,ntqxx
                do it=1,nt_max
                   if(ititpskip(it,itp)) cycle
                   ixs= ixss(it,itp) !we_ is \omega_\epsilon in Eq.(55).
                   call alagr3z2wgt(we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
                   wgt3(:,it,itp)= wfac_(it,itp)*wgt3(:,it,itp)
                   iwgt3(it,itp) = iirx(itp)*(ixs+1-2) !starting omega index ix for it,itp
                enddo                                !iirx=1 for npm=1 I think
                if(iirx(itp)/=1) call rx('sxcf: iirx=-1(TR breaking) is not yet implemented')
             enddo
             !! icount mechanism for sum ix,it,itp where W(we_(it,itp))=\sum_{i=0}^2 W(:,:,ix+i)*wgt3(i)         
             do ix = nwxi,nwx 
                do itp=lbound(zsec,1),ubound(zsec,1)
                   nit_(itp,ix)=count([((.not.ititpskip(it,itp)).and.iwgt3(it,itp)==ix,it=1,nt_max)])
                   if(nit_(itp,ix)/=0)write(6,ftox)'icou ix itp ncou/nall=',icount,ix,itp,nit_(itp,ix),ntqxx*nt_max
                enddo
             enddo
             ncoumx=maxval(nit_)
             allocate(itc(ncoumx,nwxi:nwx,ntqxx))!,itpc(ncoumx,nwxi:nwx),wgt3p(0:2,ncoumx,nwxi:nwx))
             do ix = nwxi,nwx 
                do itp=lbound(zsec,1),ubound(zsec,1)
                   iit=0
                   do it=1,nt_max
                      if((.not.ititpskip(it,itp)).and.iwgt3(it,itp)==ix) then
                         iit=iit+1
                         itc(iit,ix,itp)=it !it for given ix,itp possible iit=1,nit_(itp,ix)
                      endif
                   enddo
                enddo
             enddo
             !   ix-shifting whenr reading zw(:,:,ix)
             ikeep=99999
             do ix = nwxi,nwx  !Set wv3(:,:,0:2) is for ix,ix+1,ix+2
                if(sum(nit_(:,ix))==0) cycle
                do iw=0,2
                   if(ikeep+1==ix.and.iw<2)     then ; wv3(:,:,iw)=wv3(:,:,iw+1) 
                   elseif(ikeep+2==ix.and.iw<1) then ; wv3(:,:,iw)=wv3(:,:,iw+2)
                   else
                      read(ifrcw(kx),rec=iw+ix-nw_i+1) zw ! direct access Wc(omega) = W(omega) - v
                      ww=>zw(1:ngb,1:ngb)
                      wv3(:,:,iw)=(ww+transpose(dconjg(ww)))/2d0 !hermite part
                   endif
                enddo ! wv3 should contain zw(ix+0:ix+2) now
                ikeep=ix
                do itp=lbound(zsec,1),ubound(zsec,1)
                   do iit=1,nit_(itp,ix) !for it for given itp,ix
                      it =itc(iit,ix,itp)  !wv33 gives interpolated value of W(we_(it,itp))
                      if(it <ns1 .or. ns2<it ) cycle !xxxxxxxxxxxxxxxx
                      wv33 = wv3(:,:,0)*wgt3(0,it,itp) &
                           + wv3(:,:,1)*wgt3(1,it,itp) &
                           + wv3(:,:,2)*wgt3(2,it,itp) 
                      zmelcww(itp,it,:)= zmelcww(itp,it,:) +  matmul(zmelc(itp,it,:),wv33)
                   enddo
                enddo
             enddo
           Endblock CorrR2
         EndBlock CorrelationSelfEnergyRealAxis
         if(timemix) call timeshow(" End of CorrelationSelfEnergyRealAxis:")
         call matma(zmelcww,[(transpose(zmel(:,:,itpp)),itpp=lbound(zmel,3),ubound(zmel,3))],zsec,&
              size(zsec,1), size(zmelcww,2)*size(zmelcww,3), size(zsec,2)) !  zsec accumulated
       EndBlock zmelcww
       forall(itp=lbound(zsec,1):ubound(zsec,1)) zsec(itp,itp)=dreal(zsec(itp,itp))+img*min(dimag(zsec(itp,itp)),0d0) !enforce Imzsec<0
       call Deallocate_zmel()
       deallocate(zmelc)
       if(timemix) call timeshow("   end icount do 3030")
3030 enddo MAINicountloop
  end subroutine sxcf_scz


  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,&
       ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt, we_,wfac_,ixss,ititpskip,iirx)
    implicit none
    intent(in)::               nctot,esmr,omega,ekc,freq_r,nw_i,nw,&
         ntqxx,nt0m,nt0p,ef,nwx,nwxi,nt_max,wfaccut,wtt
    intent(out)::                                      we_,wfac_,ixss,ititpskip,iirx
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
                call rx( ' sxcf: ixs+1>nwx yyy2a')
             endif
             if(omg < ef .and. abs(ixs+1)> abs(nwxi) ) then
                call rx( ' sxcf: ixs-1<nwi yyy2b')
             endif
          endif
       enddo
    enddo
  end subroutine weightset4intreal
  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,&
       nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,&
       nwxi,nwx,nt_max)
    implicit none
    intent(in)::     omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,&
         nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug
    intent(out)::     nwxi,nwx,nt_max
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
311    enddo !continue
301 enddo !continue               
    if(nw_i==0) then          !time reversal
       nwxi = 0
       nwx  = max(ixsmx+1,ixsmin+1)
    else                      !no time revarsal working?
       nwxi = -ixsmin-1
       nwx  =  ixsmx+1
    endif
    if(nwx > nw .or. nwxi < nw_i ) call rx( ' get_nwx : |w-e| > max(w)')
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
401 enddo !continue               
  end subroutine get_nwx
end module m_sxcfsc
