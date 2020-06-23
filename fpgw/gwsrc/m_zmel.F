!> Get the matrix element zmel =  ZO^-1 <MPB psi|psi> , where ZO is ppovlz.
!! To use this module, set data in this module, and call "call get_zmelt" or "call get_zmelt2". 
!! Then we have matrix elements zmel (exchange=F for correlation)
!! or zmeltt (exchange=T). In future, they may be unified...
      module m_zmel
!! Base data for crystal structure.
!! these are set by 'call genallcf_v3' usually in the main routine.
      use m_genallcf_v3,only: 
     i     nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     i     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw,
     i     alat,delta,deltaw,esmr,symgrp,iclass,nlnmv, !,diw,dw
     i     invg, nlnmc,         !nindx,konf
     i     icore,ncore,occv,unoccv ,
     i     occc,unoccc, nocc, nunocc, plat, pos,z,ecore,  symgg,
     i     done_genallcf_v3  
!! Get eigenfuncitons. cphi is coefficients of MTO+lo part, geig is IPW parts.
!! Before calling them (get coefficients of eigen funcitons), 
!! We need to call init_readeigen, init_readeigen2 in main rouitne.
      use m_readeigen,only: Readcphi,Readgeig 
!! these are set by 'call rdpp' 
      use m_rdpp,only:       
     i     nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr,
     i     done_rdpp  
!! BZ data. To set these data 'call read_BZDATA'
      use m_read_bzdata,only: 
     i     nqbz,nqibz,  qbas,ginv,qbz,qibz,wbz,
     i     done_read_bzdata  
      use m_keyvalue,only: Getkeyvalue ! general purpose routine to read values in GWinput file.
      use m_readhbe,only: nband
      use m_itq,only: itq,ntq
      use m_readQG,only: ngpmx,ngcmx
      use m_readclasst,only: Readclasst,iclasst !iclass and iclasst is the same?
      implicit none

!!--- OUTPUT:  zmel for exchange=F, zmeltt for exchange=T.
      complex(8),allocatable,protected :: zmel(:,:,:) !,zmeltt(:,:,:) !output
      complex(8),allocatable,protected :: zmeltt(:,:,:) !output

!! Local private ---------------------
!!    These are set by mptauof in 'call mptauof_zmel'
      integer,allocatable,private :: miat(:,:)     ! mptauof_zmel
      real(8),allocatable,private :: tiat(:,:,:),shtvg(:,:)  
      real(8),allocatable,private :: ppbir(:,:,:)  ! ppbafp_v2_zmel 
      complex(8),allocatable,private :: ppovlz(:,:)        
      real(8),private:: qbasinv(3,3),q_bk(3)=1d10,qk_bk(3)=1d10
      logical,private:: init=.true.
      complex(8),allocatable,private :: cphiq(:,:), cphim(:,:)
      real(8),allocatable,private :: rmelt(:,:,:),cmelt(:,:,:)

      contains
      subroutine Dconjg_zmel()
      zmel = dconjg(zmel)
      end
      subroutine Deallocate_zmel()
      deallocate(zmel)
      end
      subroutine Deallocate_zmeltt()
      deallocate(zmeltt)
      end
!
!     set ppovlz for given q      
      subroutine setppovlz_sxcf(q,ngc,zcousq,ngb)
      intent(in)::              q,ngc,zcousq,ngb
      integer :: ngc,ngb
      real(8) :: q(3)
      complex(8) :: zcousq(ngb,ngb)
      complex(8),allocatable :: ppovl(:,:)
      if(allocated(ppovlz)) deallocate(ppovlz)
      allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
      call readppovl0(q,ngc,ppovl)
      ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
      ppovlz(nbloch+1:nbloch+ngc,:)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
      deallocate(ppovl)
      end subroutine setppovlz_sxcf
      
!! set ppovlz for given q      
      subroutine setppovlz(q,ngc,zcousq,ngb)
      intent(in)::         q,ngc,zcousq,ngb
      integer :: ngc,ngb
      real(8) :: q(3)
      complex(8) :: zcousq(ngb,ngb)
      complex(8),allocatable :: ppovl_(:,:),ppovl(:,:),ppovlzinv(:,:)
      logical:: eibz4x0
      integer:: i,nbloch
      nbloch=ngb-ngc !number of MTOs. !ngc is the number of APWs
      if(allocated(ppovlz)) deallocate(ppovlz)
      if(allocated(ppovlzinv)) deallocate(ppovlzinv)
      if(allocated(ppovl)) deallocate(ppovl)
      allocate( ppovl(ngc,ngc),ppovlz(ngb,ngb),   ppovlzinv(ngb,ngb))
      call readppovl0(q,ngc,ppovl) !q was qq
      allocate(ppovl_(ngb,ngb))
      ppovl_=0d0
      do i=1,nbloch
        ppovl_(i,i)=1d0
      enddo
      ppovl_(nbloch+1:nbloch+ngc,nbloch+1:nbloch+ngc)=ppovl
      if(.not.eibz4x0()) then   !sep2014 added for eibz4x0=F
         ppovl_= matmul(ppovl_,zcousq)
      endif
      ppovlz = ppovl_
      deallocate(ppovl_,ppovl)
      end subroutine setppovlz
      
!! set ppovlz
      subroutine setppovlz_ee(q,ngc,ngb,zcousq)
      intent(in)::            q,ngc,ngb,zcousq
      integer:: ngc,ngb
      real(8):: q(3)
      complex(8):: zcousq(:,:)
      integer:: nbloch
      complex(8),allocatable :: ppovl(:,:)
      if(allocated(ppovlz)) deallocate(ppovlz)
      if(allocated(ppovl)) deallocate(ppovl)
      allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb)) !,   ppovlzinv(ngb,ngb))
      call readppovl0(q,ngc,ppovl)
      nbloch = ngb-ngc
      ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
      ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
      end subroutine setppovlz_ee
      
!! set ppovlz
      subroutine setppovlz_chipm(zzr,ngb,nmbas1)
      intent(in)::               zzr,ngb,nmbas1
      integer::ngb,nmbas1
      complex(8):: zzr(ngb,nmbas1)
      if(allocated(ppovlz)) deallocate(ppovlz)
      allocate(ppovlz(ngb,nmbas1))
      ppovlz= zzr
      end subroutine setppovlz_chipm
      
!! Set miat,tiat,invgx,shtvg
      subroutine mptauof_zmel(symops,ng)
      intent(in)::            symops,ng
      integer:: ng
      real(8):: symops(9,ng)
      integer,allocatable ::  invgx(:)
      call Readclasst() !true class from CLASS file (genalloc_v3 may contain iclass=iatom)
      allocate(invgx(ng),miat(natom,ng),tiat(3,natom,ng),shtvg(3,ng))
      call mptauof(symops,ng,plat,natom,pos,iclasst
     o ,miat,tiat,invgx,shtvg ) !note: miat,tiat,shtvg are defined in m_zmel.
      deallocate(invgx)
      end subroutine mptauof_zmel

!! Set ppbir
      subroutine ppbafp_v2_zmel (ng,nspin,
     i                 il,in,im, mnl,
     d                 nl,nn,nclass,mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,lmxax,   !rotated CG
     i                 ppbrd)           !radial integrals
      intent(in)::               ng,nspin,
     i                 il,in,im, mnl,
     d                 nl,nn,nclass,mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,lmxax,   !rotated CG
     i                 ppbrd           !radial integrals
      integer :: ng,nspin,nclass,mnlx,mdimx,nl,nn
      integer :: il(*),in(*),im(*), lx(*),nx(*),mnl(*)
      integer :: nxx,lmxax
      real(8) :: cgr(*), ppbrd(*)
      integer:: is,irot
      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,ng,nspin))
      do irot = 1,ng
      do is = 1,nspin
         call ppbafp_v2 (irot,ng,is,nspin,
     i                 il,in,im, mnl,
     d                 nl,nn,nclass,mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,lmxax,   !rotated CG
     i                 ppbrd,           !radial integrals
     o        ppbir(:,irot,is)) !this is in m_zmel
      enddo  
      enddo  
      end subroutine ppbafp_v2_zmel

!     ! ------------------------------------
      subroutine get_zmelt(exchange,q,kx, kvec,irot,rkvec,kr,isp,  ngc,ngb,nmmax,nqmax, nctot,ncc)
      intent(in)::         exchange,q,kx, kvec,irot,rkvec,kr,isp,  ngc,ngb,nmmax,nqmax, nctot,ncc
!! Get <phiq(q,ncc+nqmax,ispq) |phim(q-rkvec,nctot+nmmax,ispm) MPB(rkvec,ngb)> ZO^-1
!!
!! ncc=0
!! kvec is in the IBZ, rk = Rot_irot(kvec), kx,kr are dummy.
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
      logical:: exchange
      integer:: kx,kr,isp,ngc,ngb,nmmax,nqmax,irot,ispq,ispm,nmini,nqini, nctot,ncc
      real(8) ::  quu(3),q(3), kvec(3),rkvec(3)
      ispq = isp
      ispm = isp
      nmini=1
      nqini=1
      call get_zmelt2(exchange,  
     &   kvec,irot,rkvec,ngc,ngb,    !MPB        for  MPB_rkvec
     &      nmini,nmmax,ispm,nctot,  !middle-phi for  phi_{q-rkvec}
     &    q,nqini,nqmax,ispq,ncc  )  !end-phi    for  phi_q
      end subroutine get_zmelt
!! ----------------------------------------
cold  ntqxx--->nqmax
cold  nbmax -->nmmax
!!note: For usual correlation mode, I think nctot=0
!!note: For self-energy mode;   we calculate <iq1|\Sigma |iq2> , where iq1 and iq2 are in nqmax.
!!       nstate = nctot+nmmax
!!       allocate(zmeltt(ngb,    nstate,  nqmax))
!!       zmeltt= < MPB     phi   | phi   > (but true matrix elements are for <phi|phi MPB> (complex conjugate).
!!               <rkvec q-rkvec  |  q    >
 !                      cphim    | cphiq 
!                       ispm     | ispq
!            nctot+  nmini:nmmax | ncc + nqini:nqmax
!                    middle state| end state
!
!!--- For dielectric funciton, we use irot=1 kvec=rkvec=q. We calulate \chi(q).
!!              q      rkvec     | q + rkvec  
!                    nkmin:nkmax | nkqmin:nkqmax
!                   (we fix nkmin=1)
!           or
!              nt0=nkmax-nkmin+1 | ntp0=nkqmax-nkqmin+1
!                      1:nt0     | 1:ntp0 
!                         occ    | unocc     
!                      (cphi_k   | cphi_kq !in x0kf)
!                    middle state| end state
!
!!     rkvec= rk(:,k)-qq  ! <phi(q+rk,nqmax)|phi(rk,nctot+nmmax)  MPB(q,ngb )>
!!     kvec = rk(:,k)-qq  ! k
!!
!! NOTE: dimension
!!  nmtot = nctot+ nmmax-mnini+1
!!  nqtot = ncc  + nqmax-nqini+1
!! <rkvec,1:ngb   q-rkvec, 1:nmtot | q, 1:nqtot>
!! ----------------------------------------
      subroutine get_zmelt2(exchange,
     &   kvec,irot,rkvec,ngc,ngb,  !MPB for       MPB_rkvec
     &   nmini,nmmax,ispm,nctot,   !middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq,ncc)  !end state for phi_q
      intent(in)::          exchange,
     &   kvec,irot,rkvec,ngc,ngb,  !MPB for       MPB_rkvec
     &   nmini,nmmax,ispm,nctot,   !middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq,ncc   !end state for phi_q
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
      logical:: exchange
      integer:: invr,nxx,itp,irot,isp,kr,no,nmmax,ngc,ngb,nqmax,nbcut
      integer:: iatomp(natom),nmini,nqini,nctot,ncc
      real(8) :: symope(3,3),shtv(3),tr(3,natom),qk(3),det
     &  , quu(3),q(3), kvec(3),rkvec(3),wtt
      complex(8),allocatable :: zzzmel(:,:,:),zw (:,:)
      integer:: nmtot,nqtot
      real(8),allocatable :: drealzzzmel(:,:,:), dimagzzzmel(:,:,:) ,ppb(:)
      logical:: debug=.false.
      complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
      complex(8):: expikt(natom)
      integer:: it,ia,kx,verbose,nstate,imdim(natom)
      logical:: oncew
      real(8),parameter::tolq=1d-8
      integer::ispq,ispm,iii,itps
!!
      if(allocated(zmel)) deallocate(zmel)
      if(allocated(zmeltt)) deallocate(zmeltt)
      if(verbose()>80) debug=.true.
      if(debug) write(*,*) 'get_zmel in m_zmel: start'
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
      if(.not.done_genallcf_v3) call rx('m_zmel: not yet call genallcf_v3')
      if(.not.done_rdpp)        call rx('m_zmel: not yet call rdpp')
      if(.not.done_read_bzdata) call rx('m_zmel: not yet call read_bzdata')

      if(init) then
         call minv33(qbas,qbasinv)
         allocate( cphiq(nlmto,nband), cphim(nlmto,nband))
         init=.false.
      endif

      if(sum(abs(q-q_bk))>tolq) then
         print *,' go to readcphi',q,nlmto,ispq,quu
ccccccccccccccccccccccccccccccccc         
         call readcphi(q, nlmto,ispq, quu, cphim )
         cphiq(1:nlmto,1:ntq) = cphim(1:nlmto,itq(1:ntq))
         q_bk=q
      endif

      allocate( rmelt(ngb, nctot+nmmax, ncc+nqmax), ! nstate= nctot+nband
     &  cmelt(ngb, nctot+nmmax, ncc+nqmax))
      if(debug) write(*,*) 'get_zmel in m_zmel: 22222222'

!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      qk =  q - rkvec        
      if(sum(abs(qk-qk_bk))>tolq) then
         call readcphi(qk, nlmto,ispm, quu, cphim)
         qk_bk=qk
      endif
c      call getsrdpp2( nclass,nl,nxx)
!! Rotate atomic positions invrot*R = R' + T
      invr  =  invg(irot)       !invrot (irot,invg,ngrp)
      tr    = tiat(:,:,invr)
      iatomp= miat(:,invr)
      symope= symgg(:,:,irot)
      shtv  = matmul(symope,shtvg(:,invr))
!! ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
!! Note spin-dependence. Look for ixx==8 in hbas.m.F calling basnfp.F, which gives ppbrd.
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass))
      ppb = ppbir(:,irot,ispq) 
      if(debug) write(*,*) 'get_zmel in m_zmel: 3333333333'

!TIME1_1001 "init"
!TIME0_1101

!! phase factors expikt(ia) is for exp(ik.T(R))
      do ia = 1,natom
        imdim(ia)  = sum(nblocha(iclass(1:ia-1)))+1
        expikt(ia) = exp(img *tpi* sum(kvec*tr(:,ia)) ) 
      end do
      nmtot  = nctot + nmmax -nmini+1      ! = phi_middle
      nqtot  = ncc   + nqmax -nqini+1      ! = phi_end
      allocate( zzzmel(nbloch,nmtot,nqtot))
      zzzmel=0d0
!! MTO Core 
      if(ncc>0.or.nctot>0) then
        call psicb_v3  ( nctot,ncc,nmmax,nqmax,iclass,expikt,
     i              cphim(1,nmini),  !middle phi
     i              cphiq(1,nqini),  !end phi
     i              ppb,!ppb,            
     i              nlnmv,nlnmc,nblocha, !mdim, 
     i              imdim,iatomp,
     i              mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i              icore,ncore,nl,nnc,
     o              zzzmel)
      endif
      if(debug) write(6,'("Goto psi2b_v3 nctot ncc nmmax nqmax=",4i4)') nctot,ncc,nmmax,nqmax
      if(debug) write(6,'("4444 zzzmelsum ",3i5,3d13.5)') nbloch,nmtot,nqtot,sum(abs(zzzmel)),sum(zzzmel)
!! MTO Valence
      if(nmmax*nqmax>0) then      ! val num of nm  ! val num of nq
        call psi2b_v3( nctot,ncc, nmmax-nmini+1,   nqmax-nqini+1, iclass,expikt, !phase,
     i              cphim(1,nmini), 
     i              cphiq(1,nqini),
     i              ppb,! ppb,         
     i              nlnmv, nlnmc,nblocha, !mdim,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx, natom,nclass,
     o              zzzmel)
      endif  
      if(debug) write(6,'("5555 zzzmelsum ",3i5,3d13.5)') nbloch,nmtot,nqtot,sum(abs(zzzmel)),sum(zzzmel)
!TIME1_1101 "psi2b_v3"

!TIME0_1201 
!! IPW 
      allocate(drealzzzmel(nbloch,nmtot,nqtot),dimagzzzmel(nbloch,nmtot,nqtot))
      drealzzzmel=dreal(zzzmel)
      dimagzzzmel=dimag(zzzmel)
      deallocate(zzzmel)
!    qk =  q - rkvec   !ncc+nqmax? nqtot?
      itps = nqini
      call drvmelp( q, nqmax-nqini+1, ! q     nt0 (in FBZ)
     i  qk,  nmmax-nmini+1,            ! q-rk  ntp0
     i  kvec,        ! k in IBZ for mixed product basis. rk = symope(kvec) 
     i  ispq,ispm,ginv,
     i  ngc,ngcmx, ngpmx,nband,itq, 
     i  symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i  drealzzzmel, dimagzzzmel, nbloch, nctot,ncc,itps,
     o  rmelt,cmelt)
      if(debug) write(6,*) ' sxcf_fal1: end of drvmelp2 sum rmelt cmelt',sum(rmelt),sum(cmelt)
      deallocate(drealzzzmel,dimagzzzmel)
      if(verbose()>50) call timeshowx("5 after drvmelp")
      if(nbcut/=0.and.(.not.exchange)) then
         do it= nctot+1,nctot+min(nbcut,nmmax)
            rmelt(:, it,:) =0d0
            cmelt(:, it,:) =0d0
         enddo
      endif
!TIME1_1201 "drvmelp"

!! NOTE:=========================================
!! zmelt = rmelt(igb(rkvec), iocc(q), iunocc(q-rkvec)) + i* cmelt
!! iunocc: band index at target  q.
!! iocc:   band index at intermediate vector qk = q - rkvec 
!! igb: index of mixed product basis       at rkvec (or written as rk)
!!   igb=1,ngb 
!!   ngb=nbloch+ngc  ngb: # of mixed product basis
!!                   nbloch: # of product basis (within MTs)
!!                   ngc: # of IPW for the Screened Coulomb interaction.
!!                   igc is for given 
!! See readgeig in drvmelp2.
!! ==================================================
c-------------------------------------------------------------------------
c$$$!! smbasis 
c$$$!! smbasis ---need to fix this
!!  Read pomatr
c$$$      if(smbasis()) then  !this smbasis if block is from hsfp0.sc.m.F
c$$$         write(6,*)' smooth mixed basis : augmented zmel'
c$$$         ifpomat = iopen('POmat',0,-1,0) !oct2005
c$$$         nkpo = nqibz+nq0i
c$$$         nnmx=0
c$$$         nomx=0
c$$$         do ikpo=1,nkpo
c$$$            read(ifpomat) q_r,nn_,no,iqx !readin reduction matrix pomat
c$$$            if(nn_>nnmx) nnmx=nn_
c$$$            if(no>nomx) nomx=no
c$$$            allocate( pomat(nn_,no) )
c$$$            read(ifpomat) pomat
c$$$            deallocate(pomat)
c$$$         enddo
c$$$         isx = iclose("POmat")
c$$$         ifpomat = iopen('POmat',0,-1,0) !oct2005
c$$$         allocate( pomatr(nnmx,nomx,nkpo),qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
c$$$         do ikpo=1,nkpo
c$$$            read(ifpomat) qrr(:,ikpo),nn_,no,iqx !readin reduction matrix pomat
c$$$            nnr(ikpo)=nn_
c$$$            nor(ikpo)=no
c$$$            read(ifpomat) pomatr(1:nn_,1:no,ikpo)
c$$$         enddo
c$$$         isx = iclose("POmat")
c$$$         write(6,*)"Read end of POmat ---"
c$$$      endif
c-------------------------------------
c$$$               if(smbasis()) then !
c$$$                  ntp0= nqmax
c$$$                  nn= nnr(kx)
c$$$                  no= nor(kx)
c$$$                  allocate( pomat(nn,no) )
c$$$                  pomat= pomatr(1:nn,1:no,kx)
c$$$                  if( sum(abs(kvec-qrr(:,kx)))>1d-10 .and.kx <= nqibz ) then
c$$$                     call rx( 'qibz/= qrr')
c$$$                  endif
c$$$                  if(no /= ngb.and.kx <= nqibz) then
c$$$!!     A bit sloppy check only for kx<nqibz because qibze is not supplied...
c$$$                     write(6,"(' q  ngb  ',3d13.5,3i5)")  kvec,ngb
c$$$                     write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
c$$$                     call rx( 'x0kf_v2h: POmat err no/=ngb')
c$$$                  endif
c$$$                  if(timemix) call timeshow("xxx2222 k-cycle")
c$$$                  ngb = nn      ! Renew ngb !!!
c$$$                  allocate ( zmel  (nn, nctot+nmmax, ntp0) )
c$$$                  call matm( pomat, dcmplx(rmelt,cmelt), zmel, 
c$$$     &                 nn, no, (nctot+nmmax)*ntp0 )
c$$$                  deallocate(rmelt, cmelt)
c$$$                  allocate( rmelt(ngb, nctot+nmmax, ntp0), !ngb is reduced.
c$$$     &                 cmelt(ngb, nctot+nmmax, ntp0) )
c$$$                  rmelt = dreal(zmel)
c$$$                  cmelt = dimag(zmel)
c$$$                  deallocate(zmel,pomat)
c$$$               else
c$$$                  nn=ngb
c$$$                  no=ngb
c$$$               endif

c      if( oncew() ) then
c         write(6,"('ngb nn no=',3i6)") ngb,nn,no
c      endif
c               if(timemix) call timeshow("22222 k-cycle")
      if(allocated(zzzmel))deallocate(zzzmel) !rmel,cmel)
      if(debug) write(6,*) ' sxcf: goto wtt'
      if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)
!! === End of zmelt ; we now have matrix element zmelt= rmelt + img* cmelt === 
!TIME0_1301

!! Multipled by ppovlz and reformat
      if(exchange) then
         if(debug) write(*,*) 'exchange mode 0000 ngb nmtot nqtot',ngb,nmtot,nqtot
         allocate( zmel (ngb, nmtot, nqtot))
         zmel = dcmplx (rmelt,cmelt)
         if(debug) write(*,*) 'exchange mode 1111'
         deallocate(rmelt,cmelt)
         if(debug) then
            do it = 1,nmtot
               write(6,"('wwwwwsc ',i5,2f10.4)") it,sum(abs(zmel(:,it,1)))
            enddo
            write(*,*) 'eeeeeeeeeeeee end of wwwwsc',nctot,nmmax
            write(6,*)'sumcheck ppovlz=',sum(abs(ppovlz(:,:)))
         endif
!! OUTPUT zmeltt for exchange
         allocate(zmeltt(nmtot,nqtot,ngb))

         if(verbose()>39) then
            write(*,*)'info: USE GEMM FOR SUM (zmeltt=zmel*ppovlz) in sxcf_fal2.sc.F'
            write(*,*)'zgemmsize',nqtot*nmtot,ngb,ngb
            write(*,*)'size ,zmel',size(zmel,dim=1),size(zmel,dim=2),size(zmel,dim=3)
            write(*,*)'size ,ppovlz',size(ppovlz,dim=1),size(ppovlz,dim=2)
            write(*,*)'size ,zmeltt',size(zmeltt,dim=1),size(zmeltt,dim=2),size(zmeltt,dim=3)
         endif
         call flush(6)	
         call zgemm('T','N',nqtot*nmtot,ngb,ngb,(1d0,0d0),
     .     zmel,ngb,ppovlz,ngb,(0d0,0d0),zmeltt,nqtot*nmtot )
         deallocate(zmel)
      else
!! Correlation case. Get zmel
         if(debug) write(*,*) 'correlation mode 0000'
c         nstate = nctot + nmmax ! = nstate for the case of correlation
         allocate(zmeltt(ngb, nmtot, nqtot))
         zmeltt= dcmplx (rmelt,-cmelt) !zmeltt= <itp|it,ib>
         deallocate(rmelt,cmelt)
!! zmel(igb,it*itp) = C(ppovlz)*N(zmeltt(:,it*itp))
!!  C means Hermitian conjugate, N means normal
!! http://www.netlib.org/lapack/lapack-3.1.1/html/zgemm.f.html
!! OUTPUT
         allocate( zmel (ngb, nmtot, nqtot) )

         if(debug) write(6,'("4 zzzppp222aaa ",3d13.5)') sum(abs(zmeltt)),sum(zmeltt)
         call zgemm('C','N',ngb, nmtot*nqtot,ngb,(1d0,0d0),
     .     ppovlz, ngb, zmeltt,ngb, (0d0,0d0),zmel,ngb)
         deallocate(zmeltt)
         if(debug) write(*,*)'zz000 nmtot,ngb,nstate ',nmtot,ngb,nqtot
         if(debug) write(*,*)'zz000 sumchk zmel ',sum(abs(zmel(1:ngb,1:nmtot,1:nqtot)))
         if(debug) write(*,*) 'correlation mode end'
!TIME1_1301 "matmul_zmelp_povlz"
      endif         
      end subroutine get_zmelt2
      end module m_zmel

      
      subroutine timeshowx(info)
      character*(*) :: info
      write(*,'(a,$)')info
      call cputid(0)
      end

