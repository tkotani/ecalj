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
module m_sxcf_main
  use m_readqg,only   : Readqg0
  use m_readeigen,only: Readeval
  use m_keyvalue,only   : Getkeyvalue
  use m_zmel,only : Get_zmel_init, Setppovlz, Setppovlz_chipm, Deallocate_zmel, zmel
  use m_itq,only: itq,ntq,nbandmx
  use m_genallcf_v3,only: nlmto,nspin,nctot,niw,ecore !,symgg
  use m_read_bzdata,only: qibz,qbz,wk=>wbz,nqibz,nqbz,wklm,lxklm,nq0i, wqt=>wt,q0i, irk
  use m_readVcoud,only:   Readvcoud, vcoud,ngb,ngc
  use m_readfreq_r,only: freq_r, nw_i,nw,freqx,wx=>wwx,nblochpmx,mrecl,expa_,npm,nprecx
  use m_rdpp,only: Rdpp, nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
  use m_readqg,only:  ngpmx,ngcmx
  use m_readhbe,only: nband,mrecg
  use m_readgwinput,only: ua_,  corehole,wcorehole
  use m_mpi,only: MPI__sxcf_rankdivider
  use m_ftox
  use m_sxcf_count,only: ncount,ispc,kxc,irotc,ipc,krc,nstateMax,nstti,nstte,nstte2, nwxic, nwxc
!  use m_hamindex,only: ngrp
!  use m_eibzhs,only: nrkip=>nrkip_all,irkip_all
!  use m_read_bzdata,only: qbz,nqibz,ginv,irk,nqbz
  implicit none
  public sxcf_scz_main, zsecall
  complex(8),allocatable,target:: zsecall(:,:,:,:) !output
  private
contains
  subroutine sxcf_scz_main(ef,esmr,exchange,ixc,nspinmx) 
    intent(in)             ef,esmr,exchange,ixc,nspinmx
    logical :: exchange
    integer :: isp,nspinmx,jobsw  !nqibz, nbandmx(nqibz,nspinmx)
    real(8) :: ef,esmr !, qvec(3,nqibz)
    real(8):: ebmx
    complex(8),pointer::zsec(:,:), ww(:,:)
    integer,allocatable :: ifrcw(:),ifrcwi(:)
    integer :: ip, it, itp, i, ix, kx, irot, kr
    integer :: nt0p, nt0m,nstate , nbmax, ntqxx 
    integer :: nt,iw,ivc,ifvcoud,ngb0
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
    logical,parameter :: debug=.false.,timemix=.false.
    logical ::   oncew, onceww !, eibz4sig  
    real(8),allocatable:: we_(:,:),wfac_(:,:)
    complex(8),allocatable:: w3p(:),wtff(:)
    logical:: tote=.false.!, hermitianW
!    real(8),allocatable:: vcoud_(:),wfft(:)
    logical:: iprx,cmdopt0
    integer:: ixx,ixc,icount,ns1,ns2,ns1c,ns2c,ns1v
    real(8),parameter:: pi=4d0*datan(1d0), fpi=4d0*pi, tpi=8d0*datan(1d0),ddw=10d0
    integer:: kxold,nccc,icount0
    complex(8),allocatable:: zmelc(:,:,:)
    integer,allocatable::ndiv(:),nstatei(:,:),nstatee(:,:),irkip(:,:,:,:)
    integer:: ns2r
    iqini = 1
    iqend = nqibz             
    allocate(zsecall(ntq,ntq,nqibz,nspinmx)) 
    zsecall = 0d0
    if(.not.exchange) then! Read WV* containing W-v in MPB
       allocate(ifrcw(iqini:iqend),ifrcwi(iqini:iqend))
       do kx=iqini,iqend
          if(any(kx==kxc)) then !only for requied files for this rank
             open(newunit=ifrcw(kx), file='WVR.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
             open(newunit=ifrcwi(kx),file='WVI.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
          endif
       enddo
    endif
    kxold=-9999     !nccc= max(ncount/50,1) !just for printing
    MAINicountloop: do 3030 icount=1,ncount   !we only consider bzcase()==1 now
       !write(6,*)'do 3030 icount=',icount
       isp =ispc(icount)
       kx  =kxc(icount) !for W(kx), kx is irreducible !kx is main axis
       irot=irotc(icount)
       ip  =ipc(icount)
       kr  =krc(icount) !=irkip(isp,kx,irot,ip) runs all the k mesh points required for q(ip).
       nwxi=nwxic(icount) !minimum omega for W
       nwx =nwxc(icount)  !max omega for W
       qibz_k = qibz(:,kx)
       q(1:3)= qibz(1:3,ip)
       eq = readeval(q,isp) !readin eigenvalue
       omega(1:ntq) = eq(1:ntq)  !1:ntq
       qbz_kr= qbz (:,kr)     !rotated qbz vector. 
       qk =  q - qbz_kr       !<M(qbz_kr) phi(q-qbz_kr)|phi(q)>
       ekq = readeval(qk, isp) 
       ekc(1:nctot)= ecore(1:nctot,isp) ! core
       ekc(nctot+1:nctot+nband) = ekq (1:nband)
       nt0  = count(ekc<ef) 
       nt0p = count(ekq<ef+ddw*esmr) +nctot 
       nt0m = count(ekq<ef-ddw*esmr) +nctot
       ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
       zsec => zsecall(1:ntqxx,1:ntqxx,ip,isp)
       wtt = wk(kr)
       !if(eibz4sig()) wtt=wtt*nrkip(isp,kx,irot,ip)
       ns1 =nstti(icount) ! Range of middle states is [ns1:ns2] 
       ns2 =nstte(icount) ! 
       ns2r=nstte2(icount) !Range of middle states [ns1:ns2r] for CorrelationSelfEnergyRealAxis
!       write(6,ftox)'do3030:isp kx ip',isp,kx,ip,'icou/ncou=',icount,ncount,'ns1:ns2=',ns1,ns2
       ZmelBlock:block !zmel= <M(qbz_kr,ib) phi(it,q-qbz_kr,isp) |phi(itpp,q,isp)> 
       !                       ib=1,ngb !MPB,  it=ns1:ns2 ! MiddleState, itpp=1:ntqxx ! EndState, 
         if(kxold/=kx) then
            call Readvcoud(qibz_k,kx,NoVcou=.false.) !Readin ngc,ngb,vcoud ! Coulomb matrix
            call Setppovlz(qibz_k,matz=.true.)       !Set ppovlz overlap matrix used in Get_zmel_init in m_zmel
            if(debug) write(6,*) ' sxcf_fal2sc: ngb ngc nbloch=',ngb,ngc,nbloch
            kxold =kx
         endif
         call Get_zmel_init(q,qibz_k,irot,qbz_kr, ns1,ns2,isp, 1,ntqxx,isp, nctot,ncc=0,iprx=debug)
                            !Return zmel middle=> ns1:ns2
       endblock ZmelBlock
       Exchangemode: if(exchange) then      
          ExchangeSelfEnergy: Block
            real(8):: wfacx,vcoud_(ngb),wtff(ns1:ns2) !range of middle states ns1:ns2
            vcoud_= vcoud
            if(kx == iqini) vcoud_(1)=wklm(1)*fpi*sqrt(fpi)/wk(kx) !voud_(1) is effective v(q=0) in the Gamma cell.
            wtff= [(wfacx(-1d99, ef, ekc(it), esmr),it=ns1,ns2)]
            do concurrent(itp=lbound(zsec,1):ubound(zsec,1),itpp=lbound(zsec,2):ubound(zsec,2))
               block
                 complex(8):: w3p(ns1:ns2)
                 w3p(ns1:ns2)=[(sum(dconjg(zmel(:,it,itp))*vcoud_(:)*zmel(:,it,itpp)),it=ns1,ns2)]
                 ns1v= max(nctot+1,ns1) ! minimum index of valence for core+valence
                 w3p(ns1v:ns2) = w3p(ns1v:ns2) * wtff(ns1v:ns2)
                 if(corehole) then !not checked well
                    ns2c= min(ns2,nctot) !ns2c upper limit of core index of core+valence
                    ns1c= ns1 !ns1c lower limit of core index of core+valence
                    w3p(ns1c:ns2c) = w3p(ns1c:ns2c) * wcorehole(ns1c:ns2c,isp)
                 endif
                 zsec(itp,itpp) = zsec(itp,itpp) - wtt * sum( w3p(:) )
               end block
            enddo
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
       zmelcww: Block !range of middle states is [ns1:ns2] 
         complex(8):: zmelcww(1:ntqxx,ns1:ns2,1:ngb)
         allocate(zmelc(1:ntqxx,ns1:ns2,1:ngb)) 
         forall(itp=1:ntqxx) zmelc(itp,:,:)=transpose(dconjg(zmel(:,:,itp))) 
         CorrelationSelfEnergyImagAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
           use m_purewintz,only: wintzsg_npm_wgtim
           real(8):: esmrx(ns1:ns2), wgtim(0:npm*niw,ntqxx,ns1:ns2)
           complex(8),target:: zw(nblochpmx,nblochpmx),zwz(ns1:ns2,ntqxx,ntqxx)
           logical:: init=.true.
           if(timemix) call timeshow(" CorrelationSelfEnergyImagAxis:")
           esmrx(ns1:nctot)= 0d0
           esmrx(max(nctot+1,ns1):ns2) = esmr
           itpdo:do concurrent(itp = lbound(zsec,1):ubound(zsec,1), it=ns1:ns2)
              call wintzsg_npm_wgtim(npm,ua_,expa_,we=.5d0*(omega(itp)-ekc(it)),esmr=esmrx(it),  &
                   wgtim=wgtim(:,itp,it))! Integration weight wgtim along im axis for zwz(0:niw*npm)
           enddo itpdo
           zmelcww=0d0
           iwimag:do ixx=0,niw !concurrent may(or maynot) be problematic because zmelcww is for reduction 
!           iwimag:do concurrent(ixx=0:niw) !niw is ~10. ixx=0 is for omega=0 nw_i=0 (Time reversal) or nw_i =-nw
              if(ixx==0) then ! at omega=0 ! nw_i=0 (Time reversal) or nw_i =-nw
                 read(ifrcw(kx),rec=1+(0-nw_i)) zw ! direct access read Wc(0) = W(0) - v
              elseif(ixx>0) then ! 
                 read(ifrcwi(kx),rec=ixx) zw ! direct access read Wc(i*omega)=W(i*omega)-v
              endif
              ww=>zw(1:ngb,1:ngb)
              associate( nleft=>size(zmelc,1)*size(zmelc,2), nm=>size(ww,1), nright=>size(ww,2))
                call matmaw(zmelc,ww,zmelcww, nleft, nm, nright, wtt*wgtim(ixx,:,:)) !pure
              end associate
           enddo iwimag
         EndBlock CorrelationSelfEnergyImagAxis
         CorrelationSelfEnergyRealAxis: Block !Fig.1 PHYSICAL REVIEW B 76, 165106(2007)
           real(8):: we_(ns1:ns2r,ntqxx),wfac_(ns1:ns2r,ntqxx)
           integer:: ixss(ns1:ns2r,ntqxx),iirx(ntqxx)
!           logical:: ititpskip(ns1:ns2r,ntqxx)
           complex(8),target:: zw(nblochpmx,nblochpmx)
           if(timemix) call timeshow(" CorrelationSelfEnergyRealAxis:")
           if(debug)write(6,*)' CorrelationSelfEnergyRealAxis: Block '
           call weightset4intreal(nctot,esmr,omega,ekc,freq_r,nw_i,nw,&
                ntqxx,nt0m,nt0p,ef,nwx,nwxi,ns1,ns2r,wfaccut,wtt,&
                we_,wfac_,ixss,iirx)
           if(any(iirx(1:ntqxx)/=1)) call rx('sxcf: iirx=-1(TR breaking) is not yet implemented')
           CorrR2:Block
             real(8):: wgt3(0:2,ns1:ns2r,ntqxx),amat(3,3)!3-point interpolation weight for we_(it,itp) 
             complex(8)::zadd(ntqxx),wv33(ngb,ngb),wv3(ngb,ngb,0:2)
             integer:: iwgt3(ns1:ns2r,ntqxx),i1,i2,iw,ikeep,ix
             integer:: nit_(ntqxx,nwxi:nwx),icountp,ncoumx,iit,irs
             integer,allocatable:: itc(:,:,:),itpc(:,:)
             if(timemix) call timeshow(" CorrR2:")
             do concurrent( itp=1:ntqxx, it=ns1:ns2r) !it=ns1:ns2) !itp:end states, it:middle states
                !we_ is \omega_\epsilon in Eq.(55).
                associate(ixs => ixss(it,itp) , x=>we_(it,itp),xi=>freq_r(ixs-1:ixs+1)) 
                  if(ixs==0) cycle
                  !call alagr3z2wgt(we_(it,itp),freq_r(ixs-1),wgt3(:,it,itp))
                  amat(1:3,1) = 1d0
                  amat(1:3,2) = xi(1:3)**2
                  amat(1:3,3) = xi(1:3)**4
                  wgt3(:,it,itp)= wfac_(it,itp)*matmul([1d0,x**2,x**4], inverse33(amat)) 
                  iwgt3(it,itp) = iirx(itp)*(ixs+1-2) !starting omega index ix for it,itp
                endassociate
             enddo
             ! icount mechanism for sum ix,it,itp where W(we_(it,itp))=\sum_{i=0}^2 W(:,:,ix+i)*wgt3(i)         
             do concurrent(ix = nwxi:nwx,itp=lbound(zsec,1):ubound(zsec,1))
                nit_(itp,ix)=count([(ixss(it,itp)>0.and.iwgt3(it,itp)==ix, it=ns1,ns2r)])
             enddo
             ncoumx=maxval(nit_) 
             allocate(itc(ncoumx,nwxi:nwx,ntqxx)) !,itpc(ncoumx,nwxi:nwx),wgt3p(0:2,ncoumx,nwxi:nwx))
             do concurrent(ix=nwxi:nwx, itp=1:ntqxx) !lbound(zsec,1):ubound(zsec,1))
                block
                  integer:: iit,it
                  iit=0
                  do it=ns1,ns2r
                     if(ixss(it,itp)>0.and.iwgt3(it,itp)==ix) then
                        iit=iit+1
                        itc(iit,ix,itp)=it !it for given ix,itp possible iit=1,nit_(itp,ix)
                     endif
                  enddo
                  nit_(itp,ix)=iit
                endblock
             enddo
             !   ix-shifting whenr reading zw(:,:,ix)
             ikeep=99999
             do ix = nwxi,nwx  !Set wv3(:,:,0:2) is for ix,ix+1,ix+2
                if(sum(nit_(:,ix))==0) cycle
                if(ikeep+1==ix) then ! use wv3 at previous ix. a cash mechanism
                   wv3(:,:,0)=wv3(:,:,1) 
                   wv3(:,:,1)=wv3(:,:,2)
                   irs=2
                elseif(ikeep+2==ix) then 
                   wv3(:,:,0)=wv3(:,:,2) 
                   irs=1
                else !ikeep+n==ix where n>2
                   irs=0
                endif   
                do iw=irs,2
                   read(ifrcw(kx),rec=iw+ix-nw_i+1) zw ! direct access Wc(omega) = W(omega) - v
                   ww=>zw(1:ngb,1:ngb)
                   wv3(:,:,iw)=(ww+transpose(dconjg(ww)))/2d0 !hermite part
                enddo ! wv3 should contain zw(ix+0:ix+2) now
                ikeep=ix
                do itp=lbound(zsec,1),ubound(zsec,1)
                   do iit=1,nit_(itp,ix) !for it for given itp,ix
                      it =itc(iit,ix,itp)  !wv33 gives interpolated value of W(we_(it,itp))
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
         block
           integer:: n1,nm,n2
           complex(8):: zmeltr(ns1:ns2,1:ngb,lbound(zsec,2):ubound(zsec,2))
           do concurrent (itpp=lbound(zsec,2):ubound(zsec,2))
              zmeltr(:,:,itpp)=transpose(zmel(:,:,itpp))
           enddo
           n1=size(zsec,1); nm=size(zmelcww,2)*size(zmelcww,3); n2=size(zsec,2)
           zsec= zsec+ matmul(reshape(zmelcww,[n1,nm]),reshape(zmeltr,[nm,n2]))
         endblock
       EndBlock zmelcww
       forall(itp=lbound(zsec,1):ubound(zsec,1)) zsec(itp,itp)=dreal(zsec(itp,itp))+img*min(dimag(zsec(itp,itp)),0d0) !enforce Imzsec<0
       call Deallocate_zmel()
       deallocate(zmelc)
       if(timemix) call timeshow("   end icount do 3030")
3030 enddo MAINicountloop
  endsubroutine sxcf_scz_main
  subroutine weightset4intreal(& ! generate required data set for main part of real part integration.
         nctot,esmr,omega,ekc,freq_r,nw_i,nw,ntqxx,nt0m,nt0p,ef,nwx,nwxi,ns1,ns2r,wfaccut,wtt,&
         we_,wfac_,ixss,iirx)
    implicit none
    intent(in)::&
         nctot,esmr,omega,ekc,freq_r,nw_i,nw,ntqxx,nt0m,nt0p,ef,nwx,nwxi,ns1,ns2r,wfaccut,wtt
    intent(out)::&
         we_,wfac_,ixss,iirx !,ititpskip
    integer:: ntqxx,nctot,nw_i,nw,nt0m,nwx,nwxi,ns2r,ns1
    real(8):: ef,omega(ntqxx),ekc(ntqxx),freq_r(nw_i:nw),esmr,wfaccut,wtt
    real(8):: we_(ns1:ns2r,ntqxx),wfac_(ns1:ns2r,ntqxx)
    integer:: ixss(ns1:ns2r,ntqxx),iirx(ntqxx)
!    logical:: ititpskip(ns1:ns2r,ntqxx)
    integer:: itini,iii,it,itend,wp,itp,iwp,nt0p,ixs
    real(8):: omg,esmrx,wfacx2,we,wfac,weavx2
    ixss=0
    do itp = 1,ntqxx          !this loop should finish in a second
       omg = omega(itp)
       iirx(itp) = 1
       if( omg < ef .and. nw_i/=0) iirx(itp) = -1
       if (omg >= ef) then
          itini= max(ns1,nt0m+1)
          itend= ns2r
          iii=  1
       else
          itini= max(1,ns1)
          itend= min(nt0p,ns2r)
          iii= -1
       endif
       do it = itini,itend     ! nt0p corresponds to efp
!          ititpskip(it,itp)=.false.
          esmrx = esmr
          if(it<=nctot) esmrx = 0d0
          wfac_(it,itp) = wfacx2(omg,ef, ekc(it),esmrx)
          wfac = wfac_(it,itp)
          if(wfac<wfaccut) then
!             ititpskip(it,itp)=.true.
             cycle 
          endif
          wfac_(it,itp)=  wfac_(it,itp)*wtt*iii
          !   Gaussian smearing we_= \bar{\omega_\epsilon} in sentences next to Eq.58 in PRB76,165106 (2007)
          !   wfac_ = $w$ weight (smeared thus truncated by ef). See the sentences.
          we_(it,itp) = .5d0* abs( omg-weavx2(omg,ef, ekc(it),esmr) ) 
          we= we_(it,itp) 
          do iwp = 1,nw 
             ixs = iwp
             if(freq_r(iwp)>we) exit
          enddo
          ixss(it,itp) = ixs
       enddo
    enddo
  end subroutine weightset4intreal
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
  pure function inverse33(matrix) result(inverse) !Inverts 3X3 matrix
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
