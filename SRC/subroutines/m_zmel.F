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
     &         Setppovlz, Setppovlz_chipm,  Mptauof_zmel ,GramSchmidt_zmel,rwzmel, drvmelp3,ppbafp_v2
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
      subroutine rwzmel(iq,k,isp_k,isp_kq,rw)
      intent(in)::      iq,k,isp_k,isp_kq,rw
      integer:: iq,k,isp_k,isp_kq,izmel,nmtot_,nqtot_
      character(1):: rw
      character*8 :: charext
      open(newunit=izmel,file='zmel.'//trim(charext(iq))//'_'
     &     //trim(charext(k))//'_'//trim(charext(isp_k))//trim(charext(isp_kq)),
     &     form='unformatted')
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
c$$$  end subroutine

      
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine ppbafp_v2 (ig,ng,isp,nspin,
     i                 il,in,im, mnl,
     d                 nl,nn,nclass,mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,lmxax,   !rotated CG
     i                 ppbrd,           !radial integrals
     o  ppb)
c calculates for all R (class)
c <Phi(RLn) Phi(RL'n') B(R,i)>
c   n differenciate core phi phidot localOrbital.
c
c B(R,i)   = orthonormal basis, i=1,mdim
c
c---------------------
c iln            = s. indxln.f
c nn,nl          = number of n and l
c
c mdim           = dimension of orthonormal basis
c mnl   =  (m,n,l) composite idex --- numbering.
c in,il,im      = index for n,l,m s. indxlnm.f
c
c ppb            = <Phi(RLn) Phi(RL'n') B(R,i)>
!! this rouitne show index system (ibas,l,m,n)
      implicit none
      integer(4),intent(in) :: ig,ng,isp,nspin,nclass,mnlx,mdimx
      integer(4),intent(in) :: il(mnlx,nclass),in(mnlx,nclass),im(mnlx,nclass)
      integer(4),intent(in) :: lx(nclass),nx(0: 2*(nl-1),nclass)
      integer(4),intent(in) :: nxx,lmxax
      real(8), intent(out) :: ppb(mnlx,mnlx,mdimx,nclass)
      real(8), intent(in) :: cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng)
      real(8), intent(in) :: ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nxx,nclass*nspin)
      integer(4) :: ic,nl,nn, i,lb,nb,mb,lmb,i1,ibas,i2
      integer(4) :: np,lp,mp,lmp,n,l,m,lm, mnl(nclass)
      integer(4):: zvztest
      do ic  = 1, nclass
        ibas = ic
!! loop over product basis
        i = 0 !i = product basis index.
        do lb  = 0, lx (ibas)
          do nb  = 1, nx (lb,ibas)
            do mb  = -lb, lb
              i    = i+1  !The number of product basis is  =(i at the end of loop).
              lmb  = lb*lb + lb + mb + 1
!! loop over LMTO basis
              do  i2 = 1,mnl(ic)
                np   = in(i2,ic)
                lp   = il(i2,ic)
                mp   = im(i2,ic)
                lmp  = lp*lp + lp + mp + 1

                do  i1 = 1,mnl(ic)
                  n    = in(i1,ic)
                  l    = il(i1,ic)
                  m    = im(i1,ic)
                  lm   = l*l + l + m + 1
                  ppb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
     &                  * ppbrd( l,n, lp,np, lb,nb, isp+nspin*(ic-1))
                  if(zvztest()==2) then
                    if(abs(ppb(i1,i2,i,ic))/=0d0  )then
                      write(6,"('i1 i2 i ic=',4i4,' ppb=',f10.6)") i1,i2,i,ic, ppb(i1,i2,i,ic)
                    endif
                  endif
cccccccccccccccccccccccccccccccccccccccccccccccccc
c        pdb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(1,2))
c        dpb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(2,1))
c        ddb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(2,2))
ccccccccccccccccccc
c          write(6,'(4i4,2d15.6,7i2)')lm,lmp,i,lmb, cgr(lm,lmp, lmb,ig),
c     &       ppbrd( l,n, lp,np, lb,nb, ioff(1,1))
c     &       ,l,n, lp,np, lb,nb, ioff(1,1)
c         endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      end
c--------------------------------------------------------------------

c$$$c-------------------------------------------------------
c$$$      subroutine drvmelp( iq, ntp0, kp,nt0, k,
c$$$     i        ngc,ngcmx,ngpn, ngpmx,ngvecpB,ngveccBr,geigB,nband,itp0,
c$$$     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
c$$$     i       rmel, cmel, nbloch,noccx,nctot,
c$$$     o       rmelt,cmelt)
c$$$c ! q1=q    ---> iq 1:ntp0 q-point
c$$$c ! q2=q-rk ---> kp 1:nt0  occupied
c$$$c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c$$$c  zmelpl
c$$$c-----------------------------------------------------
c$$$      implicit none
c$$$      integer(4):: iq,kp,k,ngp1, ngp2, ngpmx,nqbz,nqibz,
c$$$     &      ngcmx ,nctot,nband,
c$$$     &      ntp0,nt0,nbloch,noccx,  itx, ngc,nnum,inum,ig1,ig2,igc,
c$$$     &             ngvecpB  (3,ngpmx,nqbz),
c$$$     &             ngveccBr (3,ngcmx,nqibz), itp0(ntp0),
c$$$     &     nadd(3),ngpn(nqbz)
c$$$      complex(8),allocatable::  zmelpl(:,:,:),geigq(:,:)
c$$$      real(8):: qbas(3,3),shtv(3),qdiff(3),add(3)
c$$$     &         ,qibz(3,nqibz),qbz(3,nqbz),qbasinv(3,3),symope(3,3)
c$$$     &         ,rmel(nbloch,noccx,ntp0)
c$$$     &         ,cmel(nbloch,noccx,ntp0)
c$$$     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
c$$$     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
c$$$     & ,pi=3.1415926535897932D0
c$$$      complex(8) :: geigB  (ngpmx,nband,nqbz)
c$$$c      write(6,*)' drvmelp',ngpn(iq),ngpn(kp),ngc
c$$$      ngp1 = ngpn(iq)  ! q        in FBZ
c$$$      ngp2 = ngpn(kp)  ! q-rk     in FBZ
c$$$c     ngc  = ngcni(k)  ! k-points in IBZ
c$$$      allocate( geigq(ngp1,ntp0), zmelpl(ngc,nt0,ntp0) )
c$$$      do itx = 1,ntp0
c$$$        geigq(1:ngp1,itx) = geigB(1:ngp1,itp0(itx), iq)
c$$$      enddo
c$$$c     write(6,*)' drvmelp x1',symope,qibz(1:3,k),qbz(:,iq),qbz(:,kp)
c$$$      qdiff = matmul(symope,qibz(1:3,k)) - qbz(:,iq) + qbz(:,kp)
c$$$      ! rk          -q          +(q-rk) is not zero.
c$$$c     write(6,*)' drvmelp qdiff',qdiff
c$$$      add  = matmul(qbasinv,qdiff)
c$$$      nadd = idint( add + dsign(.5d0,add))
c$$$c      write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)')
c$$$c     &  add,nadd,abs(add-nadd)
c$$$      if(sum(abs(add-nadd))>1d-7) stop "drvmlp: abs(add-nadd))>1d-10"
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      if(.false.) then !old version
c$$$c      call melpln( ngp1, ngvecpB  (:,1:ngp1, iq)
c$$$c     &           , ngp2, ngvecpB  (:,1:ngp2, kp)
c$$$c     &           , ngc,  ngveccBr (:,1:ngc , k ), nadd,
c$$$c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c$$$c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c$$$c     i       shtv, matmul(symope,qibz(1:3,k)), qbas,
c$$$c     o       zmelpl)
c$$$c      else             !new version from Sep-10 2000
c$$$      call melpln2( ngp1, ngvecpB  (:,1:ngp1, iq)
c$$$     &           ,  ngp2, ngvecpB  (:,1:ngp2, kp)
c$$$     &           ,  ngc,  ngveccBr (:,1:ngc , k ), nadd,
c$$$     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c$$$     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c$$$     i       shtv, matmul(symope,qibz(1:3,k)),qibz(1:3,k), symope, qbas,
c$$$     o       zmelpl)
c$$$c       endif
c$$$      deallocate(geigq)
c$$$      rmelt=0d0
c$$$      cmelt=0d0
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      zmelpl=0d0
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
c$$$     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
c$$$      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
c$$$     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      write(6,*)' !!!!!! NO PLANE WAVE !!!'
c$$$c      return
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c extensiton to nbloch+ngc
c$$$      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c$$$     & = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
c$$$      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c$$$     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
c$$$      deallocate(zmelpl)
c$$$
c$$$      end

      
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
!> Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
!! q1=q    ---> iq 1:ntp0 q-point
!! q2=q-rk ---> kp 1:nt0  occupied
!!  zmelpl
      subroutine drvmelp( q, ntp0, q_rk,nt0, qik, ispq,ispq_rk,ginv,
     i        ngc,ngcmx,ngpmx_in,nband,itp0, 
     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,nctot,ncc,itps,
     o       rmelt,cmelt)
      use m_readqg, only : readqg
      use m_readeigen,only : readgeigf
c-----------------------------------------------------
      implicit none
      real(8),intent(in):: q(3),q_rk(3),qik(3),ginv(3,3)
      integer(4),intent(in):: ngpmx_in,nqbz,nqibz,ngcmx ,nctot,nband,
     &  ntp0,nt0,nbloch,ngc,itp0(ntp0),ispq,ispq_rk,itps,ncc
      real(8),intent(in):: qbas(3,3),shtv(3),qibz(3,nqibz),qbz(3,nqbz),
     &  qbasinv(3,3),symope(3,3),rmel(nbloch,nctot+nt0,ncc+ntp0),
     &  cmel(nbloch,nctot+nt0,ncc+ntp0)
      real(8),intent(out):: rmelt(1:nbloch+ngc,1:nctot+nt0,1:ncc+ntp0)
      real(8),intent(out):: cmelt(1:nbloch+ngc,1:nctot+nt0,1:ncc+ntp0)

      integer(4) :: ngp1, ngp2
      integer(4) :: ngvecpB1(3,ngpmx_in),ngvecpB2(3,ngpmx_in),nadd(3)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      real(8) :: qdiff(3),add(3),tolq=1d-5
      complex(8) :: geig1(ngpmx_in,nband),geig2(ngpmx_in,nband)
      complex(8),allocatable::  zmelpl(:,:,:)
      integer(4)::verbose
      logical:: debug=.false.
      if(verbose()>=90) debug=.true.
      if(verbose()>50) write(*,*) ' drvmelp: start'
      call readqg('QGpsi',q,    qt,   ngp1, ngvecpB1) !qt is used q in practice.
      call readqg('QGpsi',q_rk, q_rkt,ngp2, ngvecpB2)
      if(debug) write(6,*)'drvmelp: goto readgeig'
      geig1= readgeigf(q,ispq)       !call readgeig(q,    ngpmx_in, ispq, qu1, geig1)
      geig2= readgeigf(q_rk,ispq_rk) !call readgeig(q_rk, ngpmx_in, ispq_rk, qu2, geig2)
c      if(debug) write(6,*)'drvmelp2: endof readgeig'
c      if(sum(abs(qt-qu1))>1d-10) call rx( 'drvmelp2;qu1/=qu1x')
c      if(sum(abs(q_rkt-qu2))>1d-10) call rx( 'drvmelp2;qu2/=qu2x')
      if(debug) write(6,*)' drvmelp2:end of read geig '
      qdiff = matmul(symope,qik)  - qt + q_rkt
      ! rk    -q  +(q-rk) is not zero.
!! <rk q-rk |q>
      add  = matmul(qbasinv,qdiff)
      nadd = idint( add + dsign(.5d0,add))
      if(sum(abs(add-nadd))>tolq) then
        write(6,*)' drvmelp qdiff',qdiff
        write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)') 
     &  add,nadd,abs(add-nadd)
        call rx( "drvmlp: abs(add-nadd))>1d-10")
      endif
      if(debug) write(6,*)' goto alloc zmelpl'
      if(debug) write(6,*) 'drvmelp2:goto melpln2t'
      if(ngc/=0) then
        allocate( zmelpl(ngc,nt0,ntp0) )
        call melpln2t(ngp1, ngvecpB1
     &   ,  ngp2, ngvecpB2
     &   ,   ngc,  nadd,
     &   geig1(1:ngp1,itps-1+itp0(1:ntp0)), ntp0, ! q1=(shifted q) ->iq ngp1 1:ntp0 q-point
     &   geig2(1:ngp2,1:nt0), nt0, ! q2=(shifted q-rk) -> kp ngp2 1:nt0  occupied 
     i   shtv, matmul(symope,qik),qik, symope, qbas,
     i   qt,                    !qt oct2013 means q1
     o   zmelpl)
      endif
      if(debug) write(6,*)'end of melpln2t'
      rmelt=0d0
      cmelt=0d0
      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = 
     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =  
     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
!! extensiton to nbloch+ngc
      if(ngc/=0) then
        rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     &   = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
        cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     &   = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
        deallocate(zmelpl)
      endif
      if(debug) write(6,*)'end of drvmelp2'
      end


c$$$!> Get weight for each k-points near 0.
c$$$      subroutine q0iwgt(symops,ngrp,q0i,nq0i, 
c$$$     o                  wgt0)
c$$$c wgt0(irreducible-k, irotation)
c$$$      implicit none
c$$$      integer(4) :: ixx,ix,i, ngrp,ig,nq0i
c$$$      real(8)     :: q0(3,6),q0i(3,6),symops(3,3,ngrp),sym(3,3),
c$$$     &   qt(3), q0in(3,ngrp*nq0i), wgt0(nq0i,ngrp)
c$$$c      write(6,*)' q0iwgt:'
c$$$      wgt0 = 0d0
c$$$      ixx = 0
c$$$      do i = 1,nq0i
c$$$        qt = q0i(:,i)
c$$$c equivalence check
c$$$        do ig = 1,ngrp
c$$$          sym = symops(:,:,ig)
c$$$c          write(66,'(2i3,3f12.5)')i,ig,matmul(sym,qt)
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$c          do ix = 1,ixx
c$$$c            if(sum(abs(q0in(:,ix)-matmul(sym,qt)))<1d-10) then
c$$$c              goto 1111
c$$$c            endif
c$$$c          enddo
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$          ixx = ixx+1
c$$$          q0in (:,ixx) = matmul(sym,qt)
c$$$          wgt0 (i,ig)  = 1d0
c$$$c          write(6,'(2i3,3f12.5)')i,ig,q0in (:,ixx)
c$$$ 1111     continue
c$$$        enddo
c$$$      enddo
c$$$      wgt0  = wgt0 /ixx
c$$$ccccccccccccccccccccc
c$$$c     do ig=1,ngrp
c$$$c     write(6,*)' sym ig  =',ig
c$$$c     write(6,'(3f12.6)') symops(1,1:3,ig)
c$$$c     write(6,'(3f12.6)') symops(2,1:3,ig)
c$$$c     write(6,'(3f12.6)') symops(3,1:3,ig)
c$$$c      enddo
c$$$ccccccccccccccccccccc
c$$$      end

c$$$!>Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c$$$!! rmelp, cmelp
c$$$      subroutine melpln
c$$$     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc,ngvecc, nadd,
c$$$     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
c$$$     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied                        
c$$$     i       shtv,q, qbas, 
c$$$     o       zmelp)
c$$$      implicit none
c$$$      integer(4):: ngp1, ngvecp1(3,ngp1), ngp2, ngvecp2(3,ngp2),
c$$$     &     itp,it, ngc,ngvecc(3,ngc),nnum,inum,ig1,ig2,igc,ntp0,nt0,
c$$$     &     nadd(3)
c$$$      complex(8):: geigq1(ngp1,ntp0),
c$$$     &             geigq2(ngp2,nt0), zmelp(ngc,nt0,ntp0),phase(ngc),
c$$$     &             img=(0d0,1d0)
c$$$      real(8)::  q(3),qg(3),qbas(3,3),shtv(3),qdummy(3)
c$$$     & ,pi=3.1415926535897932D0
c$$$      do igc = 1,ngc
c$$$        qg(1:3) = 2*pi * (q(1:3)+ matmul(qbas, ngvecc(1:3,igc)))
c$$$        phase(igc) = exp(img*sum(qg*shtv))
c$$$      enddo
c$$$      zmelp=0d0
c$$$      do ig1=1,ngp1
c$$$        do ig2=1,ngp2
c$$$          do igc=1,ngc
c$$$            if(sum(abs( ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2)
c$$$     &             -ngvecc(1:3,igc) -nadd))  ==0) then
c$$$c          inum = inum+1
c$$$! igdelta(3,inum) means the delta function delta_G1,G2+G.
c$$$c          igdelta(1:3,inum) = (/ig1,ig2,igc/)
c$$$              do itp=1,ntp0
c$$$                do it =1,nt0
c$$$                  zmelp(igc,it,itp) = zmelp(igc,it,itp)
c$$$     &        + dconjg(geigq1(ig1,itp))*geigq2(ig2,it) * phase(igc)
c$$$                enddo
c$$$              enddo
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      end

!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
!> Mattrix elements <Plane psi |psi> from interstitial plane wave.
!! zmelp(igc(qi),it(q2),itp(q1)) = <itp(for q1+G1)| it(for q2+G2) igc> 
!! NOTE: shtv = g(delta_{g^-1})
      subroutine melpln2t
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied                        
     i       shtv,q, qi, symope, qbas, !
     i       qt, !qt oct2013 for G1
     o       zmelp)
c      use m_mpi, only:mpi__rank
      use m_read_ppovl,only: getppx2,
     &  nvggg,nvgcgp2,ngvecc, nggg,ngcgp,ngcread, ggg,ppovlinv,
     &  ngc2,ngvecc2
      implicit none
      integer(4):: itp,igc
      integer(4),intent(in) :: ngp1, ngvecp1(3,ngp1), ngp2
      integer(4),intent(in) :: ngvecp2(3,ngp2), ngc,nadd(3),ntp0,nt0
      complex(8),intent(in) :: geigq1(ngp1,ntp0),geigq2(ngp2,nt0)
      real(8),intent(in) :: shtv(3),q(3),qi(3), symope(3,3),qbas(3,3)
      real(8),intent(in) :: qt(3)
      complex(8),intent(out) :: zmelp(ngc,nt0,ntp0)
      complex(8),parameter :: img=(0d0,1d0)
      real(8), parameter :: pi=3.1415926535897932D0
      integer(4) :: nn(1:3)
      integer::igcgp2,iggg,igp1,igp2
      integer:: nxminggg,nxmaxggg,nyminggg,nymaxggg,nzminggg,nzmaxggg
      integer:: nxmingcgp2,nxmaxgcgp2,nymingcgp2,nymaxgcgp2,nzmingcgp2,nzmaxgcgp2
      integer,allocatable:: igggi(:,:,:),igcgp2i(:,:,:),ngveccR(:,:)
      complex(8),allocatable::ggitp(:,:),gp2phas2(:),phase(:)
      integer:: ngcgp2,ngcs(1)
      complex(8)::zdotc
      complex(8),allocatable:: zmelp0(:,:,:),ggitp_(:,:)
      complex(8),allocatable:: z2(:,:)
      logical:: debug=.false.
      integer:: verbose
      if(verbose()>=90) debug=.true.
      if(debug) write(6,*)' melpln2t: '
!! getppx allocate and return ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2).
!!      call getppx(qbas,  ngcx,qi,symope) ! rppovl.f

!! allocate and return foobar decleared in use m_read_ppovl,only: foobar
      if(debug) write(6,*)' melpln2t: nadd=',nadd
      call getppx2(qbas,qi) ! rppovl.f
      if(verbose()>=100) write(6,*)' melpln2t: end of getppx2'
      if(ngc/=ngcread) then
         write(6,*)qi,ngcread,ngc
         call rx( 'melpln2: ngc/= ngcx by getppx:PPOVLG')
      endif   
      ngcgp2=ngcgp
!! --- new version --- nov2013takao
!! zmelp(igc,it,itp) = <itp(for G1)|it(for G2) igc> matrix element.
!!   zmelp0(igc,it,itp) = <G1|G2 Gc'> geig^*(G1,itp) geig(G2,it) 
!!   zmelp(igc,it,itp) =   = zmelp0(Gc',it,itp) <Gc'|Gc>^-1
!!   (<Gc'|Gc>^-1 is dconjg(ppovlinv)
!! New ggg matrix <Gc |G1 G2> is introduced.
!!
!!    <Gc G2|G1> is equivalent to <-Gc+G1-G2>; described by ggg
!! Readin input
!!    ggg(1:nggg) = <Gc+G2-G1>
!!    nvggg(3,1:nggg)   for Gc+G2-G1
!!    nvgcgp2(3,ngcgp2) for Gc+G2
!!    ppovlinv(ngc,ngc) <Gc|Gc> matrix
!!
!! Old version used ppx= <Gc|Gc>^-1 <Gc|Gc2> (Gc2 span G1-G2 ) in rdata4gw_v2
!! -------------------------------------
!!  Make igggi inversion table
      nxminggg =minval(nvggg(1,1:nggg)) ;nxmaxggg=maxval(nvggg(1,1:nggg))
      nyminggg =minval(nvggg(2,1:nggg)) ;nymaxggg=maxval(nvggg(2,1:nggg))
      nzminggg =minval(nvggg(3,1:nggg)) ;nzmaxggg=maxval(nvggg(3,1:nggg))
      if(verbose()>100) write(6,"(' nggg nvggg range: ',i6, 6i4)") 
     & nggg, nxminggg,nxmaxggg,nyminggg,nymaxggg,nzminggg,nzmaxggg
      allocate( igggi(nxminggg:nxmaxggg,nyminggg:nymaxggg,nzminggg:nzmaxggg) )
      igggi = -100000
      if(verbose()>=100) write(6,*)' melpln2t: goto do igggi'
      do iggg =1,nggg
        igggi(nvggg(1,iggg),nvggg(2,iggg),nvggg(3,iggg)) = iggg
      enddo
!! ggitp(Gc+G2)= \sum_G1 <(Symope(Gc)+G2)-G1> geigq1(G1,itp)*exp(-i*G1*shtv)*exp(-i(q-Gadd)*shtv)
!! NOTE: nvgcgp2(:,igcgp2) means symope(Gc)+ G2
      allocate(ngveccR(1:3,1:ngc))
      ngcs(1) = ngc
      if(verbose()>=100) write(6,*)' melpln2t: goto rotation of ngvecc xxxx'
      call rotgvec(symope, 1, ngc, ngcs, qbas, ngvecc, 
     o    ngveccR) 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! simple case commented out. This worked (but slow) nov2013takao.
c$$$!!      write(6,*)'--- goto simple case for debug---'
c$$$      write(6,*)'fffffff: ngc ngp1 ngp2=',ngc,ngp1,ngp2
c$$$      zmelp0 = 0d0
c$$$      do igc=1,ngc
c$$$        phas =exp( img*2d0*pi*sum((q+matmul(qbas,ngveccR(:,igc)))*shtv) ) !conjg
c$$$      do igp1=1,ngp1     !for ngp1
c$$$      do igp2=1,ngp2
c$$$c        nn = ngvecp1(:,igp1)-ngvecp2(:,igp2)-ngveccR(:,igc) - nadd 
c$$$        nn = -ngvecp1(:,igp1) + ngvecp2(:,igp2) + ngveccR(:,igc) + nadd 
c$$$        iggg = igggi(nn(1),nn(2),nn(3))
c$$$        do it=1,nt0
c$$$        do itp=1,ntp0
c$$$          zmelp0(igc,it,itp) = zmelp0(igc,it,itp)
c$$$     &      +  phas*ggg(iggg)* dconjg(geigq1(igp1,itp))* geigq2(igp2,it)
c$$$        enddo
c$$$        enddo
c$$$      enddo
c$$$      enddo
c$$$      enddo
c$$$      call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
c$$$      return
ccccccccccccccccccccccccccccccccccccccccccccccccc
      if(verbose()>=100) write(6,*)' melpln2t: godo igp1 loop'
      allocate(ggitp(ntp0,ngcgp2))
      ggitp = 0d0
!$OMP parallel private(nn,iggg)
      do igp1  = 1,ngp1   !for ngp1
!$OMP do
      do igcgp2= 1,ngcgp2 !for ngc+ngp2 
         nn = ngvecp1(:,igp1)- nvgcgp2(:,igcgp2) - nadd 
         ! G1 -(Gc+G2) - Gadd !  -Gadd= -rk + qt -q_rk
c         write(6,*)' nn=',nn
         if(nn(1)<nxminggg .or. nxmaxggg<nn(1) ) cycle
         if(nn(2)<nyminggg .or. nymaxggg<nn(2) ) cycle
         if(nn(3)<nzminggg .or. nzmaxggg<nn(3) ) cycle
         iggg = igggi(nn(1),nn(2),nn(3))
         if(iggg<0) cycle
c         write(6,*)'iggg=',iggg,nggg
!! ggg(iggg) = <qt+G1 -(rk+Gc) -(q_rk+G2) >, where 
!!       gggx(1:ntp0,igcg2) = gggx(1:ntp0,igcg2)+ ggg(iggg)*geigq1_(1:ntp0,igp1)
! this is time-consuging part
         call zaxpy(ntp0,ggg(iggg), geigq1(igp1,1),ngp1, ggitp(1,igcgp2),1)
      enddo
      enddo
!$OMP end parallel 
      if(verbose()>=100) write(6,*)' melpln2t: 11111111'
      deallocate(igggi)
      ggitp = dconjg(ggitp)
!! igcp2 inverse table
      nxmingcgp2 =minval(nvgcgp2(1,1:ngcgp2)) ;nxmaxgcgp2=maxval(nvgcgp2(1,1:ngcgp2))
      nymingcgp2 =minval(nvgcgp2(2,1:ngcgp2)) ;nymaxgcgp2=maxval(nvgcgp2(2,1:ngcgp2))
      nzmingcgp2 =minval(nvgcgp2(3,1:ngcgp2)) ;nzmaxgcgp2=maxval(nvgcgp2(3,1:ngcgp2))
      allocate(igcgp2i(nxmingcgp2:nxmaxgcgp2,nymingcgp2:nymaxgcgp2,nzmingcgp2:nzmaxgcgp2))
      igcgp2i = -100000
      do igcgp2 =1,ngcgp2
        igcgp2i(nvgcgp2(1,igcgp2),nvgcgp2(2,igcgp2),nvgcgp2(3,igcgp2))=igcgp2
      enddo
      if(verbose()>=100) write(6,"(' nggg nvggg range: ',i6, 6i4)") 
     &  ngcgp2, nxmingcgp2,nxmaxgcgp2,nymingcgp2,nymaxgcgp2,nzmingcgp2,nzmaxgcgp2
!! zmelp <=  \sum_G2 ggitp(Gc+G2) geigqg2(G2))
!! note \bfr'= g (\bfr) +\delta_g  (\bfr= {\bf r})
!! mapping of function g[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
      allocate(gp2phas2(nt0),phase(ngc))
      do igc=1,ngc
        phase(igc)=exp( img*2d0*pi*sum((q+matmul(qbas,ngveccR(:,igc)))*shtv) )
      enddo  
!! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
      allocate(zmelp0(ngc,nt0,ntp0))
      zmelp0=0d0
      if(verbose()>=100) write(6,*)' melpln2t: 33333333'
!! If range of nvgcgp2 is not large enough, igcgp2 can be -100000.
c#ifdef USE_GEMM_FOR_SUM
      if(verbose()>39) write(*,*)'info: USE GEMM FOR SUM (zmelp0=geigq2*ggitp_) in ppbafp.fal.F'
!$OMP parallel private(z2,ggitp_,igcgp2,nn)
!$OMP master
!$    write(*,*)'OMP parallel ',omp_get_num_threads()
!$OMP end master
      allocate(z2(ngc,nt0))
      allocate(ggitp_(ngp2,ngcgp2))
!$OMP do 
      do itp= 1,ntp0
      do igc=1,ngc
         do igp2=1,ngp2
            nn = ngveccR(:,igc) + ngvecp2(:,igp2)
            igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
            ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
         enddo
      enddo
         call zgemm('T','N',ngc,nt0,ngp2,(1.0d0,0.0d0),ggitp_,ngp2,geigq2,ngp2,(0.0d0,0.0d0), 
     &        z2,ngc)
         zmelp0(:,:,itp)= z2
      enddo
      deallocate(z2,ggitp_)
      if(verbose()>=100) write(6,*)' melpln2t: 4444444444'
      deallocate(ngveccR,phase)
      call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
      deallocate(igcgp2i,ggitp,gp2phas2,zmelp0)
      if(verbose()>=100) write(6,*)' melpln2t: end'
      end subroutine melpln2t
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine drvmelp3( q, ntp0, q_rk,nt0, qik, isp,ginv,
     i        ngc,ngcmx,ngpmx,nband,itp0, 
     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,noccx,nctot,
     o       zmelt) 
      use m_readqg,only: readqg
      use m_readeigen,only:readgeigw
c ! q1=q    ---> iq 1:ntp0 q-point
c ! q2=q-rk ---> kp 1:nt0  occupied                        
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c  zmelpl
c-----------------------------------------------------
      implicit none
      real(8):: q(3),q_rk(3),qik(3),ginv(3,3)
      integer(4):: ngp1, ngp2, ngpmx,nqbz,nqibz,  !iq,kp,k,
     &      ngcmx ,nctot,nband,
     &      ntp0,nt0,nbloch,noccx,  itx, ngc,nnum,inum,ig1,ig2,igc,
     &             ngvecpB1(3,ngpmx),  
     &             ngvecpB2(3,ngpmx),  
     &             ngveccBr(3,ngcmx), itp0(ntp0),
     &     nadd(3),isp  !,ngpn(nqbz)
      complex(8),allocatable::  zmelpl(:,:,:),geigq(:,:)
      real(8):: qbas(3,3),shtv(3),qdiff(3),add(3)
     &         ,qibz(3,nqibz),qbz(3,nqbz),qbasinv(3,3),symope(3,3)
     &         ,rmel(nbloch,noccx,ntp0) 
     &         ,cmel(nbloch,noccx,ntp0) 
c     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
c     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     & ,pi=3.1415926535897932D0 
!      complex(8) :: geigB  (ngpmx,nband,nqbz)
      complex(8):: geig1(ngpmx,nband),geig2(ngpmx,nband)
     &            ,zmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      integer(4)::verbose
c
c      write(6,*)' drvmelp',ngpn(iq),ngpn(kp),ngc
      if(verbose()>50) write(6,*)' drvmelp2:'
c      ngp1 = ngpn(iq)  ! q        in FBZ 
c      ngp2 = ngpn(kp)  ! q-rk     in FBZ
      call readqg('QGpsi', q,    qt,   ngp1, ngvecpB1)
      call readqg('QGpsi', q_rk, q_rkt,ngp2, ngvecpB2)
      call readgeigW(q,    ngpmx, isp, qu1, geig1)
      call readgeigW(q_rk, ngpmx, isp, qu2, geig2)
c debug
c      geig1=0d0
c      geig2=0d0

      if(sum(abs(qt-qu1))>1d-10) stop 'drvmelp2;qu1/=qu1x'
      if(sum(abs(q_rkt-qu2))>1d-10) stop 'drvmelp2;qu2/=qu2x'
      if(verbose()>=100) write(6,*)' end of read geig '

c     ngc  = ngcni(k)  ! k-points in IBZ
c      allocate( geigq(ngp1,ntp0) )
c      do itx = 1,ntp0
c        geigq(1:ngp1,itx) = geigB(1:ngp1,itp0(itx), iq)
c      enddo

c     write(6,*)' drvmelp x1',symope,qibz(1:3,k),qbz(:,iq),qbz(:,kp)
c      qdiff = matmul(symope,qibz(1:3,k)) - qbz(:,iq) + qbz(:,kp) 
      qdiff = matmul(symope,qik)  - qt + q_rkt 
                           ! rk    -q  +(q-rk) is not zero.
      add  = matmul(qbasinv,qdiff)
      nadd = idint( add + dsign(.5d0,add))
      if(sum(abs(add-nadd))>1d-7) then
        write(6,*)' drvmelp qdiff',qdiff
        write(6,'(" add nadd diff =",3f12.5,3i12,3d13.6)') 
     &  add,nadd,abs(add-nadd)
        call rx( "drvmlp3: abs(add-nadd))>1d-10")
      endif
c     write(6,*)' drvmelp qdiff',qdiff
c      write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)') 
c     &  add,nadd,abs(add-nadd)

      allocate( zmelpl(ngc,nt0,ntp0) )
c      write(6,*)' zzzzzzz zmelpl',ngc,nt0,ntp0

cccccccccccccccccccccccccccccccccccccccccccccccc
c      if(.false.) then !old version
c      call melpln( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           , ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           , ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
c     i       shtv, matmul(symope,qibz(1:3,k)), qbas, 
c     o       zmelpl)
c      else             !new version from Sep-10 2000
c       call melpln2( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           ,  ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           ,  ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
c     i       shtv, matmul(symope,qibz(1:3,k)),qibz(1:3,k), symope, qbas, 
c     o       zmelpl)
c       endif
      if (nt0 .ne. ntp0) stop 'drvmelp3: nt0 .ne. ntp0'
c      write(6,*)'wmat.F: drvmelp3: goto melpln2t',nt0,itp0(1:nt0)
       call melpln2t( ngp1, ngvecpB1
     &           ,   ngp2, ngvecpB2 
c     &           ,   ngc,  ngveccBr, nadd,
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,  ! q1=q    ---> iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,itp0(1:nt0)), nt0,     ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
     i       shtv, matmul(symope,qik),qik, symope, qbas, qt, !qt is dummy...
     o       zmelpl)
c      deallocate(geig1,geig2)
c      rmelt=0d0
c      cmelt=0d0
      zmelt=0d0
c      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = 
c     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)  
c      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =  
c     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
       zmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = dcmplx
     &(rmel (1:nbloch, 1:nctot+nt0, 1:ntp0),  
     & cmel (1:nbloch, 1:nctot+nt0, 1:ntp0))
c extensiton to nbloch+ngc   
c      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c     & = dreal( zmelpl(1:ngc,         1:nt0,    
** eigenfunciton rotation PMT basis for self-ener    1:ntp0) )  
c      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )  
       zmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     &=zmelpl(1:ngc,               1:nt0,            1:ntp0) 
      deallocate(zmelpl) 
      end
      subroutine psi2b_v3(nctot,ncc,nt0,ntp0,iclass,phase,
     i                   cphik, !rbk,cbk,rhbk,chbk,
     i                   cphikq,!rbkq,cbkq,rhbkq,chbkq,
     i                   ppb, !pdb,dpb,ddb,
     i                   nlnmv,nlnmc,mdim,imdim,iatomp,
     d                   mdimx,nlmto,nbloch,nlnmx,
     d                   natom,nclass,
     o zpsi2b)
c originaly 92.03.17 by Ferdi.
c takao modified at Apr 2002(v2) Feb2006(v3).
c
c calculates <psi(k',t') | psi(k,t) B(R,i)>
c for all R
c psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
c B(R,i)   = Bloch orthonormal product basis for atom R
c psi(k,t) is stored after nctot
c
c nt0        = no. t
c ntp0       = no. t'
c coskt,sinkt= exp(ik.T)
c cphik  b(k)
c cphikq b(k')
c
c ppb        = <phi(RLn) phi(RL'n') B(R,i)>
c
c nlnmv      = number of l,n,m for valence
c nlnmc      = number of n,l,m for core states
c mdim       = number of optimal product basis functions for each class
c nctot      = total no. allowed core states
c nbloch     = total no. optimal product basis
c nlnmx      = maximum number of l,n,m
c
c zpsi2b     =  the matrix elements
C----------------------------------------------------------
      implicit none
      integer(4):: nlmto,nctot,ncc,nt0,ntp0,nlnmx,natom,nbloch,
     & ia,ic,nc,nv,nc1,ias,iap,icp,i,mdimx,nclass,itp,jp,ib
      complex(8):: cphik(nlmto,*),cphikq(nlmto,ntp0),phase(natom)
      complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc)
      integer(4),allocatable::iasx(:)
      integer(4) :: nzwork1,nzwork2,
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),
     i          imdim(natom),iatomp(natom)
      real(8)::  ppb(nlnmx,nlnmx,mdimx,nclass)
      complex(8),allocatable :: zz(:,:), zwork(:,:),zppb(:,:)
      complex(8) :: alpha,beta

c zwork is used as work array for ppb(:,:,i,ic) and for zpsi2b(ib,:,:).
      allocate( zppb(nlnmx,nlnmx),zz(nlnmx,ntp0), zwork(nt0,ntp0) ) !Work areas
      beta=0d0  ; alpha=1d0
c      allocate(rr(nlnmx,ntp),cc(nlnmx,ntp))

c check dimensions
      if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
      if( sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
      allocate(iasx(natom))
      ias = 1
      do ia = 1,natom
        iasx(ia) = ias
        ias = ias + nlnmv(iclass(ia))
      enddo
Cstop2rx 2013.08.09 kino      if(ias-1/=nlmto) stop ' psi2b_v3:sum(nlnmv)/= nlmto'
      if(ias-1/=nlmto) call rx( ' psi2b_v3:sum(nlnmv)/= nlmto')

c loop over atoms
      do  ia = 1,natom
        ic   = iclass(ia)
        nc   = nlnmc(ic)
        nv   = nlnmv(ic)
        nc1  = nc + 1
Cstop2rx 2013.08.09 kino        if (nc+ nlnmv(ic) > nlnmx) stop 'psi2b_v3: nlnmx exceeded'
        if (nc+ nlnmv(ic) > nlnmx) call rx( 'psi2b_v3: nlnmx exceeded')
c        phase= dcmplx(coskt(ia),sinkt(ia))
        ias  = iasx(ia)
        iap  = iatomp(ia)
        icp  = iclass(iap)
        do   i = 1,mdim(icp) ! loop over optimal product basis
c---------------------------------------------------
cc sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)>
cc for a given i, for all L'n' and t'
cc bkq is complex but < > is real
c1      do     itp = 1,ntp0
c1      do      jp = 1,nlnmv(ic)
c1      zz (jp,itp)=dconjg(
c1     &    sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,nc+jp,i,icp)) )
c1      end do
c1      end do
c
c2      zz(1:nv,1:ntp0) =dconjg(
c2     & matmul(  transpose(ppb(nc1:nc+nv,nc1:nc+nv,i,icp))
c2     &         ,cphikq(ias:ias+nv-1,1:ntp0)) )
c
c3        call dgemm('T','N',nv,ntp0,nv,
c3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
c3     &          dreal(cphikq(ias:ias+nv-1,1:ntp0)), nv,
c3     &          0d0,
c3     &   rr, nlnmx )
c3        call dgemm('T','N',nv,ntp0,nv,
c3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
c3     &          dimag(cphikq(ias:ias+nv-1,1:ntp0)), nv,
c3     &          0d0,
c3     &   cc, nlnmx )
          zppb(1:nv,1:nv) = ppb(nc+1:nc+nv,nc+1:nc+nv,i,icp)
          call zgemm('T','N',nv,ntp0,nv,
     &   alpha, zppb,nlnmx, cphikq(ias,1), nlmto,  beta,
     &   zz,  nlnmx )
          do itp = 1,ntp0
            do jp = 1,nv
              zz(jp,itp)= dconjg(zz(jp,itp) )
            enddo
          enddo
c----------------------------------------------------
c <psi(k+q,t') | psi(k,t) B(i)>
c1      do      it = 1,nt0
c1      do     itp = 1,ntp0
c1       zpsi2b(ib,nctot+it,itp)=
c1     &   phase * sum( zz(1:nv,itp)*cphik(ias:ias+nv-1,it) )
c1c end of t'(unoccupied)-loop
c1      end do
c1c end of t(occupied)-loop
c1      end do
c3        call zgemm('T','N',nt0,ntp0,nv,
c3     &   phase, cphik(ias:ias+nv-1,1:nt0),  nv,
c3     &          dcmplx(rr(1:nv,1:ntp0),-cc(1:nv,1:ntp0)),  nv,
c3     &          0d0,
c3    &   zpsi2b(imdim(iap)-1+i,nctot+1:nctot+nt0,1:ntp0), nt0)
          call zgemm('T','N', nt0,ntp0,nv,
     &   phase(ia), cphik(ias,1),nlmto, zz,nlnmx, beta,
     &   zwork, nt0 )
          ib = imdim(iap)-1+i
          zpsi2b(ib,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)=zwork
c------------------------------------------------------
        end do !end of optimal product basis-loop
      end do !end of atom-loop
c      deallocate(rr,cc,iasx)
      deallocate(zz,zppb,zwork,iasx)
      end


c------------------------------------------------------------------------------------
      subroutine psicb_v3 (nctot,ncc,nt0,ntp0,  iclass, phase, 
     i                   cphik,
     i                   cphikq, 
     i                   ppb,
     i                   nlnmv,nlnmc,mdim,
     i                   imdim,iatomp,
     i                   mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i                   icore,ncore,nl,nnc,  
     o            zpsi2b)
      implicit none
      intent(in)          nctot,ncc,nt0,ntp0,  iclass, phase, 
     i                   cphik,
     i                   cphikq, 
     i                   ppb,
     i                   nlnmv,nlnmc,mdim,
     i                   imdim,iatomp,
     i                   mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i                   icore,ncore,nl,nnc  
      intent(out) zpsi2b
C- Calculates <psi (k+q,t) |core(k,t) B(R,ibloch)> ---------------------
C      also   <core(k+q,t) |psi (k,t) B(R,ibloch)>
Cr B(R,i)   = Mixed basis.
Cr core(k,t)= core states
Cr <psi(k+q,t') | core(k,t) B(R,i)> = S[RLn]  cphik(RLn,k+q,t')^*  * ppb
Cr
Ci nt0        = no. k states
Ci ntp0       = no. kq states
Ci cphik coeefficients of MT part or valence eigenfunction.
Ci icore      = index for core states
Ci ncore      = no. core states in each class
Ci ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>
Ci nlnmv      = number of l,n,m for valence
Ci nlnmc      = number of l,n,m for core states
Ci mdim       = number of optimal product basis functions for each class
Ci nbloch     = total no. optimal product basis
Ci nlnmx      = maximum number of l,n,m
Co zpsi2b     =  the matrix elements
Cr coskt,sinkt= exp(ik.T)
Ch           Ferdi 92.03.17       : original for ASA
Ch           takao at Apr2002(v2) : mixed basis version 2
Ch                    Feb2006(v3) : Add case with ncc/=0
Cr
C-------------------------------------------------------------------------
      integer(4):: nl,nlmto,nbloch,nclass,natom,nctot,ncc,nt0,ntp0,mdimx,nlnmx,nnc,
     & ib,ias,ics,ia,ic,nc,nv,nc1,iap,icp,i,itp,it,icr
      integer(4)::  icore(nl*nl*nnc,nclass),ncore(nclass),
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom)
     i          ,imdim(natom),iatomp(natom)
      real(8)::   ppb(nlnmx,nlnmx,mdimx,nclass) !  coskt(natom),sinkt(natom)
      complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc),
     &             cphikq(nlmto,*), cphik(nlmto,*),phase(natom)
      integer(4):: verbose,ixx
      if(ncc/=0 .and. ncc/=nctot) then
Cstop2rx 2013.08.09 kino        stop "psicb_v3: ncc/=0 and ncc/=ncctot"
        call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
      endif
      if(sum(ncore(iclass(1:natom)))  /= nctot) 
Cstop2rx 2013.08.09 kino     &   stop "psicb_v3:sum(ncore) wrong"
     &   call rx( "psicb_v3:sum(ncore) wrong")
c      write(6,*)" psicb_v3: ncc nctot=",ncc,nctot,sum(ncore(iclass(1:natom)))
c      stop 'qqqqqqqqqqqqqq'
      zpsi2b = 0d0
c loop over atoms
      ib         = 0
      ias        = 1
      ics        = 0
      do      ia = 1,natom
        ic         = iclass(ia)
        nc         = nlnmc(ic)
        nv         = nlnmv(ic)
        nc1        = nc + 1
c       phase  =  dcmplx(coskt(ia),sinkt(ia))
c loop over optimal product basis
        iap        = iatomp(ia)
        icp        = iclass(iap)
        ib         = imdim(iap)-1
        do       i = 1,mdim(icp)
          ib         = ib + 1
          do     itp = 1,ntp0
            do      it = 1,ncore(ic)
              icr      = icore(it,ic) ! write(6,*),it,ic,icore(it,ic)
              zpsi2b(ib,ics+it,ncc+itp) = phase(ia)* 
     &    dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))
            enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop

          if(verbose()>90) write(6,*)"psicb_v3: ia i=",ia,i
          if(ncc==0) cycle

          do     itp = 1,ncore(ic)
            do      it = 1,nt0
              icr      = icore(itp,ic) ! write(6,*),it,ic,icore(it,ic)
              if(verbose()>90) write(6,*)"psicb_v3: aaa itp it ic icr=",itp,it,ic,icr
cccccccccccccccccccccccccc
              if(.false.) then
                do ixx = 0,nv-1
                  write(6,"(5i5, 5d13.6)") ixx, nc1,icr,i,icp,cphik(ias+ixx,itp)
                  write(6,"(i5,33x,5d13.6)") ixx, ppb(nc1+ixx,icr,i,icp)
                enddo
              endif
ccccccccccccccccccccccccc
              zpsi2b(ib,nctot+it,ics+itp) = dconjg(phase(ia))* 
     &    sum(cphik(ias:ias+nv-1,it)*ppb(nc1:nc+nv,icr,i,icp))
              if(verbose()>90) write(6,*)"psicb_v3: bbb itp it=",itp,it
            enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop
          if(verbose()>90) write(6,*)"psicb_v3: xxx ia i=",ia,i
        enddo !end of optimal product basis-loop
c end of atom-loop
        ias   = ias + nlnmv(ic)
        ics   = ics + ncore(ic)
      enddo
      end
      end module m_zmel
