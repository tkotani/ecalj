!! Get the matrix element zmel =  ZO^-1 <MPB psi|psi> , where ZO is ppovlz.
!!  "call get_zmel" return zmel 
!!  All dependencies (use foobar below ) are inputs (must be protected).
      module m_zmel 
      use m_genallcf_v3,only:  
     i     nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     i     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw,
     i     alat,delta,deltaw,esmr,symgrp,iclass,nlnmv, !,diw,dw
     i     invg, nlnmc,         !nindx,konf
     i     icore,ncore,occv,unoccv ,
     i     occc,unoccc, nocc, nunocc, plat, pos,z,ecore,  symgg,
     i     done_genallcf_v3,
     &     il, in, im, mnl=>nlnm , nl,nn,nlnmx
      use m_rdpp,only: Rdpp,    !"call rdpp" generate following data.
     &     nxx,lx,nx,mdimx,nbloch,cgr,ppbrd,nblocha,done_rdpp
      use m_readeigen,only: Readcphif !,Readgeigf
      use m_read_bzdata,only: 
     i     nqbz,nqibz,  qbas,ginv,qbz,qibz,wbz,
     i     done_read_bzdata  
      use m_readgwinput,only:nbcutlow_sig
      use m_readhbe,only: nband
      use m_itq,only: itq,ntq
      use m_readQG,only: ngpmx,ngcmx
      use m_readclasst,only: Readclasst,iclasst !iclass and iclasst is the same?
!! q-dependent
      use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
!!------------------------------------------------------
!! SUBROUTINES
      public:: Get_zmel_init, Get_zmel, Get_zmel_modex0,
     &         Dconjg_zmel, Deallocate_zmel, 
     &         Setppovlz, Setppovlz_chipm,  Mptauof_zmel 
!! OUTPUT:  zmel for exchange=F, zmeltt for exchange=T.
      complex(8),allocatable,protected,public :: zmel(:,:,:)   !for  Get_zmel
!!------------------------------------------------------
!! set by Mptauof_zmel in advance
      private
      integer,allocatable,private :: miat(:,:)     ! mptauof_zmel
      real(8),allocatable,private :: tiat(:,:,:),shtvg(:,:)  
      real(8),allocatable,private :: ppbir(:,:,:)  ! ppbafp_v2_zmel 
      complex(8),allocatable,private :: ppovlz(:,:)        
      real(8),private:: qbasinv(3,3),q_bk(3)=1d10,qk_bk(3)=1d10
      logical,private:: init=.true.
      complex(8),allocatable,private :: cphiq(:,:), cphim(:,:)
      real(8),allocatable,private :: rmelt(:,:,:),cmelt(:,:,:)
!! set by Get_zmel_modex0, which is called right before Get_zmel_init in x0kf_v4h
      logical:: modex0=.false. 
      integer:: nkmin, nkqmin, isp_k, isp_kq,nmtot,nqtot
      logical:: debug=.false.

      contains
!----------------------------------------------------
      subroutine Dconjg_zmel()
      zmel = dconjg(zmel)
      end
      subroutine Deallocate_zmel()
      deallocate(zmel)
      end
c      subroutine Deallocate_zmeltt()
c      deallocate(zmeltt)
c      end
!----------------------------------------------------
      subroutine setppovlz(q,matz)
      intent(in)::         q,matz
!! set ppovlz for given q      
!  Blocked matrix multiplication      
!    ppolvz(igb,ivcou)= (1    0 ) \times  zcousq(igb, ivcou)
!                       (0 ppovl)   
!    If matz=F, no multiplication by ivcou.  Thus we have ppolz(igb,igb)
c      integer :: ngc,ngb
      real(8) :: q(3)
c      complex(8) :: zcousq(ngb,ngb)
      complex(8),allocatable :: ppovl_(:,:),ppovl(:,:)!,ppovlzinv(:,:)
      logical:: eibz4x0,matz
      integer:: i!,nbloch
c      matz=.not.eibz4x0()
c     nbloch=ngb-ngc !number of PBs. !ngc is the number of IPWs
c      if(allocated(ppovlzinv)) deallocate(ppovlzinv)
      if(allocated(ppovlz)) deallocate(ppovlz)
      if(allocated(ppovl)) deallocate(ppovl)
      allocate( ppovl(ngc,ngc),ppovlz(ngb,ngb))!,   ppovlzinv(ngb,ngb))
      call readppovl0(q,ngc,ppovl) !q was qq
      if(matz) then   !sep2014 added for eibz4x0=F
         ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
         ppovlz(nbloch+1:nbloch+ngc,:)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
c         ppovlz= matmul(ppovl_,zcousq)
      else
         ppovlz=0d0
         do i=1,nbloch
            ppovlz(i,i)=1d0
         enddo
         ppovlz(nbloch+1:nbloch+ngc,nbloch+1:nbloch+ngc) = ppovl
      endif
      deallocate(ppovl)
      end subroutine setppovlz
!----------------------------------------------------
      subroutine setppovlz_chipm(zzr,nmbas1)
      intent(in)::               zzr,nmbas1
      integer::nmbas1
      complex(8):: zzr(ngb,nmbas1)
      if(allocated(ppovlz)) deallocate(ppovlz)
      allocate(ppovlz(ngb,nmbas1))
      ppovlz= zzr
      end subroutine setppovlz_chipm
!----------------------------------------------------
      subroutine mptauof_zmel(symops,ng)
!! Set miat,tiat,invgx,shtvg, and then call ppbafp_v2_zmel successively
      intent(in)::            symops,ng
      integer:: ng
      real(8):: symops(9,ng)
      integer,allocatable ::  invgx(:)
      call Readclasst() !true class from CLASS file (genalloc_v3 may contain iclass=iatom)
      allocate(invgx(ng),miat(natom,ng),tiat(3,natom,ng),shtvg(3,ng))
      call mptauof(symops,ng,plat,natom,pos,iclasst
     o ,miat,tiat,invgx,shtvg ) !note: miat,tiat,shtvg are defined in m_zmel.
      deallocate(invgx)
      call Rdpp(ng,symops)
      call ppbafp_v2_zmel(ng)
      end subroutine mptauof_zmel
!----------------------------------------------------
      subroutine ppbafp_v2_zmel(ng) !nspin,
      intent(in)::              ng 
c     i                 il,in,im, mnl,
c     d                 nl,nn,nclass,mnlx,
c     i                 mdimx,lx,nx,nxx, !Bloch wave    
c     i                 cgr,lmxax,   !rotated CG
c     i                 ppbrd)           !radial integrals
      integer :: ng 
      integer:: is,irot
      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,ng,nspin))
      do irot = 1,ng
      do is = 1,nspin
!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
!! This is general for rotated CG coefficient
!! ppbafp_v2 generates ppbir
         call ppbafp_v2 (irot,ng,is,nspin,
     i                 il,in,im, mnl,
     d                 nl,nn,nclass,nlnmx,!mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,nl-1,       !rotated CG
     i                 ppbrd,           !radial integrals
     o        ppbir(:,irot,is)) !this is in m_zmel
      enddo  
      enddo  
      end subroutine ppbafp_v2_zmel

!! ------------------------------------
      subroutine get_zmel_modex0(n1,n2,n3,n4)
      integer:: n1,n2,n3,n4
      modex0=.true.
      nkmin =n1
      nkqmin=n2
      isp_k =n3
      isp_kq=n4
      end
!! ------------------------------------
      subroutine get_zmel_init(exchange,q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc)
      intent(in)::         exchange,q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc
!! Get <phiq(q,ncc+nqmax,ispq) |phim(q-rkvec,nctot+nmmax,ispm) MPB(rkvec,ngb)> ZO^-1
!! Right after call get_zmel_init, call get_zmel or get_zmeltt to get zmel (or zmeltt).
!! kvec is in the IBZ, rk = Rot_irot(kvec)
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
      logical:: exchange
      integer:: isp,nmmax,nqmax,irot,ispq,ispm,nmini,nqini, nctot,ncc
      real(8) ::  quu(3),q(3), kvec(3),rkvec(3)
      ispq = isp
      ispm = isp
      nmini=1
      nqini=1
      if(modex0) then
         nmini= nkmin ! (k)
         nqini= nkqmin! (k)
         ispm = isp_k
         ispq = isp_kq
      endif   
      call get_zmelt2(exchange,  
     &    kvec,irot,rkvec,          ! MPB_rkvec
     &    nmini,nmmax,ispm ,nctot,  ! middle-phi for  phi_{q-rkvec}
     &    q,nqini,nqmax,ispq ,ncc  )! end-phi    for  phi_q
      end subroutine get_zmel_init
!! ----------------------------------------
cold  ntqxx--->nqmax
cold  nbmax -->nmmax
!!note: For usual correlation mode, I think nctot=0
!!note: For self-energy mode;   we calculate <iq1|\Sigma |iq2> , where iq1 and iq2 are in nqmax.
!!       nstate = nctot+nmmax
!!       allocate(zmelt(MPB,  intermediate phi nstate,  external state phi ntqxx))
!!       zmelt= < MPB     phi   | phi   > 
!!               <rkvec q-rkvec  |  q    >
 !                      cphim    | cphiq 
!                       ispm     | ispq
!            nctot+  nmini:nmmax | ncc + nqini:ntqxx
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
!! NOTE: dimension
!!  nmtot = nctot+ nmmax-mnini+1
!!  nqtot = ncc  + ntqxx-nqini+1
!! <rkvec,1:ngb   q-rkvec, 1:nmtot | q, 1:nqtot>
!! ----------------------------------------
      subroutine get_zmelt2(exchange,
     &   kvec,irot,rkvec,          !  ngcxxx,ngbxxx,  !MPB for       MPB_rkvec
     &   nmini,nmmax,ispm,nctot,   !middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq,ncc)  !end state for phi_q
      intent(in)::          exchange,
     &   kvec,irot,rkvec,          !  ngcxxx,ngbxxx,  !MPB for       MPB_rkvec
     &   nmini,nmmax,ispm,nctot,   !middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq,ncc   !end state for phi_q
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
      logical:: exchange
      integer:: invr,nxx,itp,irot,isp,kr,no,nmmax,nqmax !, ngbxxx,ngcxxx!,nbcut !ngc,ngb,
      integer:: iatomp(natom),nmini,nqini,nctot,ncc
      real(8) :: symope(3,3),shtv(3),tr(3,natom),qk(3),det
     &  , quu(3),q(3), kvec(3),rkvec(3),wtt
      complex(8),allocatable :: zzzmel(:,:,:),zw (:,:)
c      integer:: nmtot,nqtot
      real(8),allocatable :: drealzzzmel(:,:,:), dimagzzzmel(:,:,:) ,ppb(:)
      complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
      complex(8):: expikt(natom)
      integer:: it,ia,kx,verbose,nstate,imdim(natom)
      logical:: oncew
      real(8),parameter::tolq=1d-8
      integer::ispq,ispm,iii,itps
!!
      if(allocated(zmel)) deallocate(zmel)
c      if(allocated(zmeltt)) deallocate(zmeltt)
      if(verbose()>80) debug=.true.
      if(debug) write(*,*) 'get_zmel2 in m_zmel: start'
      if(.not.done_genallcf_v3) call rx('m_zmel: not yet call genallcf_v3')
      if(.not.done_rdpp)        call rx('m_zmel: not yet call rdpp')
      if(.not.done_read_bzdata) call rx('m_zmel: not yet call read_bzdata')

      if(init) then
         call minv33(qbas,qbasinv)
         allocate( cphiq(nlmto,nband), cphim(nlmto,nband))
         init=.false.
      endif

      if(sum(abs(q-q_bk))>tolq) then
c         print *,' go to readcphi',q,nlmto,ispq,quu
ccccccccccccccccccccccccccccccccc         
         cphim= readcphif(q,ispq)
         cphiq(1:nlmto,1:ntq) = cphim(1:nlmto,itq(1:ntq))
         q_bk=q
      endif

      allocate( rmelt(ngb, nctot+nmmax, ncc+nqmax), ! nstate= nctot+nband
     &  cmelt(ngb, nctot+nmmax, ncc+nqmax))
      if(debug) write(*,*) 'get_zmel2 in m_zmel: 22222222'

!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      qk =  q - rkvec        
      if(sum(abs(qk-qk_bk))>tolq) then
c         call readcphi(qk, nlmto,ispm, quu, cphim)
         cphim= readcphif(qk, ispm)
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
      if(debug) write(*,*) 'get_zmel2 in m_zmel: 3333333333'

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
c      if(verbose()>50) call timeshowx("5 after drvmelp")
      if(nbcutlow_sig/=0.and.(.not.exchange)) then
         do it= nctot+1,nctot+min(nbcutlow_sig,nmmax)
            rmelt(:, it,:) =0d0
            cmelt(:, it,:) =0d0
         enddo
      endif

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
! INPUT rmelt,cmelt,ppovlz
      block
      complex(8),allocatable :: zmeltt(:,:,:)   !for  Get_zmel
      allocate(zmeltt(ngb, nmtot, nqtot))
      zmeltt= dcmplx (rmelt,-cmelt) !zmeltt= <itp|it,ib>
      deallocate(rmelt,cmelt)
!! zmel(igb,it,itp) = C(ppovlz)*N(zmeltt(:,it*itp))
!! C means Hermitian conjugate, N means normal
      allocate( zmel (ngb, nmtot, nqtot) )
      call zgemm('C','N',ngb, nmtot*nqtot,ngb,(1d0,0d0),
     .     ppovlz, ngb, zmeltt,ngb, (0d0,0d0),zmel,ngb)
      deallocate(zmeltt)
      end block
      end subroutine
!sssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine Get_zmel()
      end subroutine
c$$$!sssssssssssssssssssssssssssssssssssssssssssssssssssss
c$$$      subroutine get_zmeltt()   !call rightafter Get_zmelt_init(exchange=T)
c$$$      implicit none
c$$$      integer::verbose
c$$$!Get zmel.  Multipled by ppovlz and reformat, INPUT rmelt,cmelt,ppovlz
c$$$!! exchnge case ordering for given rmelt,cmelt,ppovlz
c$$$      if(debug) write(*,*) 'exchange mode 0000 ngb nmtot nqtot',ngb,nmtot,nqtot
c$$$      allocate( zmel (ngb, nmtot, nqtot))
c$$$      zmel = dcmplx (rmelt,cmelt)
c$$$      if(debug) write(*,*) 'exchange mode 1111'
c$$$      deallocate(rmelt,cmelt)
c$$$!! OUTPUT zmeltt convenient for exchange
c$$$      allocate(zmeltt(nmtot,nqtot,ngb))
c$$$      if(verbose()>39) then
c$$$         write(*,*)'info: USE GEMM FOR SUM (zmeltt=zmel*ppovlz) in sxcf_fal2.sc.F'
c$$$         write(*,*)'zgemmsize',nqtot*nmtot,ngb,ngb
c$$$         write(*,*)'size ,zmel',size(zmel,dim=1),size(zmel,dim=2),size(zmel,dim=3)
c$$$         write(*,*)'size ,ppovlz',size(ppovlz,dim=1),size(ppovlz,dim=2)
c$$$         write(*,*)'size ,zmeltt',size(zmeltt,dim=1),size(zmeltt,dim=2),size(zmeltt,dim=3)
c$$$      endif
c$$$      call flush(6)	
c$$$      call zgemm('T','N',nqtot*nmtot,ngb,ngb,(1d0,0d0),
c$$$     .     zmel,ngb,ppovlz,ngb,(0d0,0d0),zmeltt,nqtot*nmtot )
c$$$      deallocate(zmel)
c$$$      end subroutine
      end module m_zmel
