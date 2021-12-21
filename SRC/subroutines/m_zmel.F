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
      public:: Get_zmel_init, Get_zmel_modex0,
     &         Dconjg_zmel, Deallocate_zmel, 
     &         Setppovlz, Setppovlz_chipm,  Mptauof_zmel ,GramSchmidt_zmel,rwzmel
!! OUTPUT:  zmel for exchange=F, zmeltt for exchange=T.
      complex(8),allocatable,protected,public :: zmel(:,:,:)   !for  Get_zmel
!!------------------------------------------------------
!! set by Mptauof_zmel in advance
      private
      integer,allocatable,private :: miat(:,:)     
      real(8),allocatable,private :: tiat(:,:,:),shtvg(:,:)  
      real(8),allocatable,private :: ppbir(:,:,:)  
      complex(8),allocatable,private :: ppovlz(:,:)        
      real(8),private:: qbasinv(3,3),q_bk(3)=1d10,qk_bk(3)=1d10
      logical,private:: init=.true.
      complex(8),allocatable,private :: cphiq(:,:), cphim(:,:),cphitemp(:,:)
      real(8),allocatable,private :: rmelt(:,:,:),cmelt(:,:,:)
!! set by Get_zmel_modex0, which is called right before Get_zmel_init in x0kf_v4h
      logical:: modex0=.false. 
      integer:: nkmin, nkqmin, isp_k, isp_kq,nmtot,nqtot,ispq_bk,ispm_bk
      logical:: debug=.false.

      contains
!!!!!!!!!!!!!!!!!!!!!!      
      subroutine GramSchmidt_zmel() !oct15 2021
!!  allocate( zmel(ngb,     nmtot,    nqtot) )
!     !                  q,   k occ it, q+k unocc itp
      integer:: igb=1,it,itt
      complex(8):: ov(nmtot),vec(nqtot),dnorm2(nmtot)
      real(8):: dnorm
c$$$      do it = 1,nmtot        ! occ
c$$$         vec(:)= zmel(igb,it,:)
c$$$         do itt= 30,nqtot
c$$$            vec(itt)=0d0
c$$$         enddo
c$$$         zmel(igb,it,:) = vec
c$$$      enddo
      do it = 1,nmtot        ! occ
         vec(:)= zmel(igb,it,:) 
         do itt = 1,it-1
            ov(itt) = sum( dconjg(zmel(igb,itt,:))*vec(:))/dnorm2(itt)
         enddo
         vec = vec - matmul(ov(1:it-1),zmel(igb,1:it-1,:))
         dnorm2(it) = sum(dconjg(vec)*vec)
         zmel(igb,it,:) = vec
cccccccccccccccccccccc
c         zmel(igb,it,nmtot+1:) = 0d0
cccccccccccccccccccccc

      enddo
      end
      
!!!!!!!!!!!!!!!!!!!!!!      
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
      subroutine rwzmel(iq,k,isp_k,rw)
      intent(in)::        iq,k,isp_k
      integer:: iq,k,isp_k,izmel,nmtot_,nqtot_
      character(1):: rw
      character*8 :: charext
      open(newunit=izmel,file='zmel.'//trim(charext(iq))//'_'
     &     //trim(charext(k))//'_'//trim(charext(isp_k)),form='unformatted')
      if(rw=='r') then
         read(izmel) nmtot_,nqtot_
         if(allocated(zmel)) deallocate(zmel)
         allocate(zmel(1:ngb,1:nmtot_,1:nqtot_))
         read(izmel) zmel
      else
         write(izmel) nmtot,nqtot
         write(izmel) zmel
      endif
      close(izmel)
      end
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
      subroutine get_zmel_init(exchange,q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc,iprx)
      intent(in)::             exchange,q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc,iprx
!! Get <phiq(q,ncc+nqmax,ispq) |phim(q-rkvec,nctot+nmmax,ispm) MPB(rkvec,ngb)> ZO^-1
!! Right after call get_zmel_init, we get zmel defined in this module.
!! kvec is in the IBZ, rk = Rot_irot(kvec)
!! \parameter all inputs
!     ! \parameter matrix <MPB psi|psi>
      logical:: exchange
      logical,optional:: iprx
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
      if(present(iprx).and.iprx) then
         print *,'ppp modex0',modex0
         write(6,"(a,l,10i4)")'pppp3:',exchange,irot,nmini,nmmax,ispm,nctot,nqini,nqmax,ispq,ncc
         write(6,"(a,10(3f9.3,x))") 'pppp3:',q, kvec,rkvec
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
!!   nmtot = nctot+ nmmax-mnini+1
!!   nqtot = ncc  + ntqxx-nqini+1
!!   <q 1:ngb,      q-rkvec, 1:nmtot | rkvec, 1:nqtot>
!!   <end state,       middle state  |  MPB          >
!!   rkvec =mutmul(symops(:,:,irot),kvec)
!! ----------------------------------------
      subroutine get_zmelt2(exchange,
     &   kvec,irot,rkvec,          ! rkvec MPB  ngcxxx,ngbxxx,  MPB_rkvec
     &   nmini,nmmax,ispm,nctot,   ! q-rkvec middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq,ncc)  ! end state for phi_q
      intent(in)::          exchange,
     &   kvec,irot,rkvec,          
     &   nmini,nmmax,ispm,nctot,   
     &    q,nqini,nqmax,ispq,ncc   
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
      logical:: exchange
      integer:: invr,nxx,itp,irot,isp,kr,no,nmmax,nqmax !, ngbxxx,ngcxxx!,nbcut !ngc,ngb,
      integer:: iatomp(natom),nmini,nqini,nctot,ncc
      real(8) :: symope(3,3),shtv(3),tr(3,natom),qk(3),det
     &  , quu(3),q(3), kvec(3),rkvec(3),wtt
      complex(8),allocatable :: zzzmel(:,:,:),zw (:,:)
      real(8),allocatable :: drealzzzmel(:,:,:), dimagzzzmel(:,:,:) ,ppb(:)
      complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
      complex(8):: expikt(natom)
      integer:: it,ia,kx,verbose,nstate,imdim(natom)
      logical:: oncew
      real(8),parameter::tolq=1d-8
      integer::ispq,ispm,iii,itps
      
      if(allocated(zmel)) deallocate(zmel)
      if(debug) write(*,*) 'get_zmel2 in m_zmel: start'
      if(.not.done_genallcf_v3) call rx('m_zmel: not yet call genallcf_v3')
      if(.not.done_rdpp)        call rx('m_zmel: not yet call rdpp')
      if(.not.done_read_bzdata) call rx('m_zmel: not yet call read_bzdata')

      if(init) then
         call minv33(qbas,qbasinv)
         allocate( cphiq(nlmto,nband), cphim(nlmto,nband), cphitemp(nlmto,nband))
         init=.false.
      endif
!! q
cc      if(sum(abs(q-q_bk))>tolq .or. ispq/=ispq_bk)  then
      if(sum(abs(q-q_bk))>tolq .or. ispq/=ispq_bk)  then
         cphitemp= readcphif(q,ispq)
         cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq))
         q_bk=q
         ispq_bk=ispq
      endif

      allocate( rmelt(ngb, nctot+nmmax, ncc+nqmax), ! nstate= nctot+nband
     &  cmelt(ngb, nctot+nmmax, ncc+nqmax))
      if(debug) write(*,*) 'get_zmel2 in m_zmel: 22222222'

!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      qk =  q - rkvec        
      if(sum(abs(qk-qk_bk))>tolq.or.ispm/=ispm_bk) then
         cphim= readcphif(qk, ispm) !all readin but we need only bands nmini:nmmax
         qk_bk= qk
         ispm_bk= ispm
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
      if(debug) write(6,'("5555 zzzmelsum ",3i5,3d13.5)')nbloch,nmtot,nqtot,sum(abs(zzzmel)),sum(zzzmel)
      if(debug) write(6,'("6666 chim      ",4i5,4d13.5)')nmini,nmmax,nqini,nqmax,sum(cphim(:,nmini:nmmax)),sum(cphiq(:,nqini:nqmax))
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
      if(allocated(zzzmel))deallocate(zzzmel) !rmel,cmel)
      if(debug) write(6,*) ' sxcf: goto wtt'
      if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)
! INPUT rmelt,cmelt,ppovlz
      block
      complex(8),allocatable :: zmeltt(:,:,:) !for  Get_zmel
      allocate(zmeltt(ngb, nmtot, nqtot))
      zmeltt= dcmplx (rmelt,-cmelt) !zmeltt= <itp|it,ib>
      deallocate(rmelt,cmelt)
!     ! zmel(igb,it,itp) = C(ppovlz)*N(zmeltt(:,it*itp))
!     ! C means Hermitian conjugate, N means normal
      allocate( zmel (ngb, nmtot, nqtot) )
      call zgemm('C','N',ngb, nmtot*nqtot,ngb,(1d0,0d0),
     .     ppovlz, ngb, zmeltt,ngb, (0d0,0d0),zmel,ngb)
      deallocate(zmeltt)
      end block
      end subroutine
!sssssssssssssssssssssssssssssssssssssssssssssssssssss
c      subroutine Get_zmel()
c      end subroutine
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
