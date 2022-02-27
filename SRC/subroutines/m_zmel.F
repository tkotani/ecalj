!! Get the matrix element zmel =  ZO^-1 <MPB psi|psi> , where ZO is ppovlz.
!!  "call get_zmel" return zmel 
!!  All dependencies (use foobar below ) are inputs (must be protected).
      module m_zmel 
      use m_genallcf_v3,only:  
     i     nclass,natom,nspin,nl,nn,nnv,nnc, 
     i     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw,
     i     alat,delta,deltaw,esmr,iclass,nlnmv, !,diw,dw
     i     nlnmc,         !nindx,konf
     i     icore,ncore,occv,unoccv ,
     i     occc,unoccc, nocc, nunocc, plat, pos,z,ecore,  
     i     done_genallcf_v3,
     &     il, in, im, mnl=>nlnm , nl,nn,nlnmx
      use m_hamindex,only: ngrp, symgg=>symops,invg=>invgx
      use m_rdpp,only: Rdpp,    !"call rdpp" generate following data.
     &     nxx,lx,nx,mdimx,nbloch,cgr,ppbrd,nblocha,done_rdpp
      use m_readeigen,only: Readcphif !,Readgeigf
      use m_read_bzdata,only: 
     i     nqbz,nqibz,  qlat,ginv,qbz,qibz,wbz,
     i     done_read_bzdata  
      use m_readhbe,only: nband
      use m_itq,only: itq,ntq
      use m_readQG,only: ngpmx,ngcmx,Readqg
      use m_hamindex0,only: Readhamindex0,iclasst 
      use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
!!------------------------------------------------------
      public:: Get_zmel_init, Get_zmel_modex0,  Dconjg_zmel, Deallocate_zmel, 
     &         Setppovlz, Setppovlz_chipm,  Mptauof_zmel ,rwzmel, !GramSchmidt_zmel (test code)
     &   drvmelp3,ppbafp_v2 !these are for old-fashioned codes in wannier/ 
!! OUTPUT:  zmel for exchange=F, zmeltt for exchange=T.
      complex(8),allocatable,protected,public :: zmel(:,:,:),zmel0(:,:,:) !for  Get_zmel
      integer:: nbb !1st dimension of zmel
!!------------------------------------------------------
!! set by Mptauof_zmel in advance
      private
      integer,allocatable,private :: miat(:,:)     
      real(8),allocatable,private :: tiat(:,:,:),shtvg(:,:)  
      real(8),allocatable,private :: ppbir(:,:,:)  
      complex(8),allocatable,private :: ppovlz(:,:)        
      real(8),private:: qlatinv(3,3),q_bk(3)=1d10,qk_bk(3)=1d10
      logical,private:: init=.true.
      complex(8),allocatable,private :: cphiq(:,:), cphim(:,:),cphitemp(:,:)
      logical:: modex0=.false. 
      integer:: nkmin, nkqmin, isp_k, isp_kq,nmtot,nqtot,ispq_bk,ispm_bk
      logical:: debug=.false.
      integer:: nbbx=0
      contains
!sssssssssssssssssssssssssssssssssss
      subroutine GramSchmidt_zmel() !oct15 2021 !experimenal
      integer:: igb=1,it,itt
      complex(8):: ov(nmtot),vec(nqtot),dnorm2(nmtot)
      real(8):: dnorm
      do it = 1,nmtot        ! occ
         vec(:)= zmel(igb,it,:) 
         do itt = 1,it-1
            ov(itt) = sum( dconjg(zmel(igb,itt,:))*vec(:))/dnorm2(itt)
         enddo
         vec = vec - matmul(ov(1:it-1),zmel(igb,1:it-1,:))
         dnorm2(it) = sum(dconjg(vec)*vec)
         zmel(igb,it,:) = vec
      enddo
      end
!sssssssssssssssssssssssssssssssssss
      subroutine Dconjg_zmel()
      zmel = dconjg(zmel)
      end
!sssssssssssssssssssssssssssssssssss
      subroutine Deallocate_zmel()
      deallocate(zmel)
      if(allocated(zmel0)) deallocate(zmel0)
      end
!sssssssssssssssssssssssssssssssssss
      subroutine setppovlz(q,matz)
      intent(in)::         q,matz
!! set ppovlz for given q      
!    ppolvz(igb,ivcou)= (1    0 ) \times  zcousq(igb, ivcou)
!                       (0 ppovl)   
!    If matz=F, no multiplication by ivcou.  Thus we have ppolz(igb,igb)
      real(8) :: q(3)
      complex(8),allocatable :: ppovl_(:,:),ppovl(:,:)!,ppovlzinv(:,:)
      logical:: eibz4x0,matz
      integer:: i
      if(allocated(ppovlz)) deallocate(ppovlz)
      if(allocated(ppovl)) deallocate(ppovl)
      allocate( ppovl(ngc,ngc),ppovlz(ngb,ngb))!,   ppovlzinv(ngb,ngb))
      call readppovl0(q,ngc,ppovl) !q was qq
      if(matz) then   !sep2014 added for eibz4x0=F
         ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
         ppovlz(nbloch+1:nbloch+ngc,:)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
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
      nbbx=nmbas1
      end subroutine setppovlz_chipm
!----------------------------------------------------
      subroutine mptauof_zmel(symops,ng)
!! Set miat,tiat,invgx,shtvg, and then call ppbafp_v2_zmel successively
      intent(in)::            symops,ng
      integer:: ng
      real(8):: symops(9,ng)
      integer,allocatable ::  invgx(:)
      call readhamindex0()
      allocate(invgx(ng),miat(natom,ng),tiat(3,natom,ng),shtvg(3,ng))
      call mptauof(symops,ng,plat,natom,pos,iclasst
     o ,miat,tiat,invgx,shtvg ) !note: miat,tiat,shtvg are defined in m_zmel.
      deallocate(invgx)
      call Rdpp(ng,symops)
      call ppbafp_v2_zmel(ng)
      end subroutine mptauof_zmel
!----------------------------------------------------
      subroutine ppbafp_v2_zmel(ng)
      intent(in)::              ng 
      integer :: ng,is,irot
      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,ng,nspin))
      do irot = 1,ng
      do is = 1,nspin
!! ppbir is rotated <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)> by rotated cg coefficients cgr
         call ppbafp_v2 (irot,ng,is,!nspin, 
c     i                 il,in,im, mnl,
c     d                 nl,nn,nclass,nlnmx,!mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,nl-1,       !rotated CG
     i                 ppbrd,           !radial integrals
     o        ppbir(:,irot,is)) !this is in m_zmel
      enddo  
      enddo  
      end subroutine ppbafp_v2_zmel
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine ppbafp_v2 (ig,ng,isp,!nspin,
c     i                 il,in,im, mnl,
c     d                 nl,nn,nclass,nlnmx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,lmxax,   !rotated CG
     i                 ppbrd,           !radial integrals
     o  ppb)
! <Phi(RLn) Phi(RL'n') B(R,i)>. Atomic part of MPB
!  n differenciates core phi phidot localOrbital.
!  in,il,im      = index for n,l,m s. indxlnm.f
!   ppb            = <Phi(RLn) Phi(RL'n') B(R,i)>
      implicit none
      integer(4),intent(in) :: ig,ng,isp!,nspin,nclass,nlnmx,mdimx
c      integer(4),intent(in) :: il(nlnmx,nclass),in(nlnmx,nclass),im(nlnmx,nclass)
      integer(4),intent(in) :: lx(nclass),nx(0: 2*(nl-1),nclass)
      integer(4),intent(in) :: lmxax,nxx,mdimx
      real(8), intent(out) :: ppb(nlnmx,nlnmx,mdimx,nclass)
      real(8), intent(in) :: cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng)
      real(8), intent(in) :: ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nxx,nclass*nspin)
      integer(4) :: ic, i,lb,nb,mb,lmb,i1,ibas,i2 !nl,nn,
      integer(4) :: np,lp,mp,lmp,n,l,m,lm!, mnl(nclass)
      do ic  = 1, nclass
        ibas = ic
        i = 0 !i = product basis index.
        do lb  = 0, lx (ibas)
        do nb  = 1, nx (lb,ibas)
        do mb  = -lb, lb
           i    = i+1           !The number of product basis is  =(i at the end of loop).
           lmb  = lb*lb + lb + mb + 1
           do  i2 = 1,mnl(ic) !phi2 index
              np   = in(i2,ic)
              lp   = il(i2,ic)
              mp   = im(i2,ic)
              lmp  = lp*lp + lp + mp + 1
              do  i1 = 1,mnl(ic) !phi1 index
                 n    = in(i1,ic)
                 l    = il(i1,ic)
                 m    = im(i1,ic)
                 lm   = l*l + l + m + 1
                 ppb(i1,i2,i,ic)=cgr(lm,lmp,lmb,ig)*ppbrd(l,n,lp,np,lb,nb,isp+nspin*(ic-1))
              enddo
           enddo
        enddo
        enddo
        enddo
      enddo
      end
!! ------------------------------------
      subroutine rwzmel(iq,k,isp_k,isp_kq,rw,izmel0)
      intent(in)::      iq,k,isp_k,isp_kq,rw,izmel0
      logical,optional::izmel0
      integer:: iq,k,isp_k,isp_kq,izmel,nmtot_,nqtot_
      character(1):: rw
      character*8 :: charext
      open(newunit=izmel,file='zmel.'//trim(charext(iq))//'_'
     &     //trim(charext(k))//'_'//trim(charext(isp_k))//trim(charext(isp_kq)),
     &     form='unformatted')
      if(rw=='r'.and.present(izmel0)) then
         print *,'fff=',trim('zmel.'//trim(charext(iq))//'_'
     &     //trim(charext(k))//'_'//trim(charext(isp_k))//trim(charext(isp_kq)))
         read(izmel) nbb,nmtot_,nqtot_
         if(allocated(zmel0)) deallocate(zmel0)
         allocate(zmel0(1:nbb,1:nmtot_,1:nqtot_))
         read(izmel) zmel0
      elseif(rw=='r') then
         read(izmel) nbb,nmtot_,nqtot_
         if(allocated(zmel)) deallocate(zmel)
         allocate(zmel(1:nbb,1:nmtot_,1:nqtot_))
         read(izmel) zmel
      else
         write(izmel) nbb,nmtot,nqtot
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
      subroutine Get_zmel_init(q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc,iprx)
      intent(in)::             q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc,iprx
!! Get zmel= <phiq(q,ncc+nqmax,ispq) |phim(q-rkvec,nctot+nmmax,ispm) MPB(rkvec,ngb)> ZO^-1
!! kvec is in the IBZ, rk = Rot_irot(kvec)
!! \parameter all inputs
!! \parameter matrix <MPB psi|psi>
!! memo =========================================
!! igb: index of mixed product basis       at rkvec (or written as rk)
!!   igb=1,ngb 
!!   ngb=nbloch+ngc  ngb: # of mixed product basis
!!                   nbloch: # of product basis (within MTs)
!!                   ngc: # of IPW for the Screened Coulomb interaction.
!!                   igc is for given 
! zmelt= <itp|it,ib>
! zmel(igb,it,itp) = transpose(dcojng(ppovlz))*zmelt(:,it*itp)
c   matmul(symgg(:,:,irot),kvec)-rkvec can have difference of reciprocal vectors.
      logical:: iprx
      integer:: isp,nmmax,nqmax,irot,ispq,ispm,nmini,nqini, nctot,ncc
      real(8) ::  quu(3),q(3), kvec(3),rkvec(3)
      real(8),parameter::tolq=1d-8
      if(allocated(zmel)) deallocate(zmel)
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
      
      ZmeltMain: Block
      use m_readeigen,only : readgeigf
      integer:: it,ia,kx,verbose,nstate,imdim(natom),nt0,ntp0,invr, iatomp(natom),ispq_rk
      real(8):: symope(3,3),shtv(3),tr(3,natom),qk(3)
      real(8),allocatable :: ppb(:)
      complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
      complex(8):: expikt(natom)
      complex(8),allocatable::  zmelt(:,:,:)
      integer(4) :: ngp1, ngp2
      integer(4) :: ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
      real(8):: q_rkt(3),qt(3),q_rk(3)!,qik(3)
      real(8) :: qdiff(3)
      complex(8) :: geig1(ngpmx,nband),geig2(ngpmx,nband)
      logical:: debug=.false.

      if(init) then
         call minv33(qlat,qlatinv)
         allocate( cphiq(nlmto,nband), cphim(nlmto,nband), cphitemp(nlmto,nband))
         init=.false.
      endif
     
      if(sum(abs(q-q_bk))>tolq .or. ispq/=ispq_bk)  then
         cphitemp= readcphif(q,ispq)
         cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq))
         q_bk=q
         ispq_bk=ispq
      endif
!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      qk =  q - rkvec        
      if(sum(abs(qk-qk_bk))>tolq.or.ispm/=ispm_bk) then
         cphim= readcphif(qk, ispm) !all readin but we need only bands nmini:nmmax
         qk_bk= qk
         ispm_bk= ispm
      endif
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
      do ia = 1,natom
        imdim(ia)  = sum(nblocha(iclass(1:ia-1)))+1
        expikt(ia) = exp(img *tpi* sum(kvec*tr(:,ia)) )  !phase expikt(ia) is for exp(ik.T(R))
      end do
      nt0= nmmax-nmini+1
      ntp0=nqmax-nqini+1
      nmtot  = nctot + nt0     ! = phi_middle
      nqtot  = ncc   + ntp0    ! = phi_end
      allocate(zmelt(1:nbloch+ngc,nmtot,nqtot))
      zmelt=0d0
!! MTO Core part
      if(ncc>0.or.nctot>0) then
        call psicb_v3  ( nctot,ncc,nmmax,nqmax,iclass,expikt,
     i              cphim(1,nmini),  !middle phi
     i              cphiq(1,nqini),  !end phi
     i              ppb,!ppb,            
     i              nblocha, !mdim, 
     i              imdim,iatomp,
     i              mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i              icore,ncore,nl,nnc,
     o              zmelt(1:nbloch,:,:))
      endif
!! MTO Valence part
      if(nmmax*nqmax>0) then      ! val num of nm  ! val num of nq
        call psi2b_v3( nctot,ncc, nmmax-nmini+1,   nqmax-nqini+1, iclass,expikt, !phase,
     i              cphim(1,nmini), 
     i              cphiq(1,nqini),
     i              ppb,! ppb,         
     i              nblocha, !mdim,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx, natom,nclass,
     o              zmelt(1:nbloch,:,:))
      endif  
!! IPW part
      q_rk=qk
      ispq_rk=ispm
      call readqg('QGpsi',q,    qt,   ngp1, ngvecpB1) !q is mapped to qt in BZ
      call readqg('QGpsi',q_rk, q_rkt,ngp2, ngvecpB2)
      geig1= readgeigf(q,ispq)       !call readgeig(q,    ngpmx_in, ispq, qu1, geig1)
      geig2= readgeigf(q_rk,ispq_rk) 
      qdiff = matmul(symope,kvec)  - qt + q_rkt ! rk -q +(q-rk) is not zero. <rk q-rk |q>
      nadd = nint(matmul(qlatinv,qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
      call melpln2t(ngp1, ngvecpB1
     &   ,  ngp2, ngvecpB2,   ngc,  nadd,
     &   geig1(1:ngp1,nqini-1+itq(1:ntp0)), ntp0, ! q1=(shifted q) ->iq ngp1 1:ntp0 q-point
     &   geig2(1:ngp2,1:nt0), nt0, ! q2=(shifted q-rk) -> kp ngp2 1:nt0  occupied 
     i   shtv, matmul(symope,kvec),kvec, symope, qlat, 
     i   qt,                    
     o     zmelt(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0))
      nbb=ngb
      if(nbbx/=0) nbb=nbbx
      allocate( zmel(nbb,nmtot, nqtot))
      call matm(dconjg(transpose(ppovlz)), dconjg(zmelt), zmel,nbb,ngb,nmtot*nqtot)
      deallocate(zmelt)
      EndBlock ZmeltMain
      end subroutine get_zmel_init

!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
!> Mattrix elements <Plane psi |psi> from interstitial plane wave.
!! zmelp(igc(qi),it(q2),itp(q1)) = <itp(for q1+G1)| it(for q2+G2) igc> 
!! NOTE: shtv = g(delta_{g^-1})
!!--- zmelp(igc,it,itp) = <itp(for G1)|it(for G2) igc> matrix element.
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
!! ggitp(Gc+G2)= \sum_G1 <(Symope(Gc)+G2)-G1> geigq1(G1,itp)*exp(-i*G1*shtv)*exp(-i(q-Gadd)*shtv)
!! NOTE: nvgcgp2(:,igcgp2) means symope(Gc)+ G2
      subroutine melpln2t
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied                        
     i       shtv,q, qi, symope, qlat, !
     i       qt, !qt oct2013 for G1
     o       zmelp)
      use m_read_ppovl,only: getppx2,igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,
     &  nvggg,nvgcgp2,ngvecc, nggg,ngcgp,ngcread, ggg,ppovlinv, ngc2,ngvecc2
      implicit none
      integer(4):: itp,igc
      integer(4),intent(in) :: ngp1, ngvecp1(3,ngp1), ngp2
      integer(4),intent(in) :: ngvecp2(3,ngp2), ngc,nadd(3),ntp0,nt0
      complex(8),intent(in) :: geigq1(ngp1,ntp0),geigq2(ngp2,nt0)
      real(8),intent(in) :: shtv(3),q(3),qi(3), symope(3,3),qlat(3,3)
      real(8),intent(in) :: qt(3)
      complex(8),intent(out) :: zmelp(ngc,nt0,ntp0)
      complex(8),parameter :: img=(0d0,1d0)
      real(8), parameter :: pi=3.1415926535897932D0
      integer(4) :: nn(1:3)
      integer::igcgp2,iggg,igp1,igp2
      integer,allocatable:: ngveccR(:,:)
      complex(8),allocatable::ggitp(:,:),gp2phas2(:),phase(:)
      integer:: ngcgp2,ngcs(1)
      complex(8)::zdotc
      complex(8),allocatable:: zmelp0(:,:,:),ggitp_(:,:)
      complex(8),allocatable:: z2(:,:)
      logical:: debug=.false.
      integer:: verbose!,nxi,nxe,nyi,nye,nzi,nze
      if(ngc==0) return
      if(verbose()>=90) debug=.true.
      call getppx2(qlat,qi) ! read and allocate PPOVL*
      if(ngc/=ngcread)call rx( 'melpln2: ngc/= ngcx by getppx:PPOVLG')!  write(6,*)qi,ngcread,ngc
      allocate(ngveccR(1:3,1:ngc))
      ngcs(1) = ngc
      call rotgvec(symope, 1, ngc, ngcs, qlat, ngvecc, 
     o    ngveccR) 
      allocate(ggitp(ntp0,ngcgp))
      ggitp = 0d0
      do igp1  = 1,ngp1   !for ngp1
      do igcgp2= 1,ngcgp !for ngc+ngp2 
         nn = ngvecp1(:,igp1)- nvgcgp2(:,igcgp2) - nadd ! G1 -(Gc+G2) - Gadd ! -Gadd= -rk + qt -q_rk
         if(nn(1)<nxi.or.nxe<nn(1).or.nn(2)<nyi.or.nye<nn(2).or.nn(3)<nzi.or.nze<nn(3)) cycle
         iggg = igggi(nn(1),nn(2),nn(3)) !inversion table
         if(iggg<0) cycle ! ggg(iggg) = <qt+G1 -(rk+Gc) -(q_rk+G2) >, where 
         ggitp(:,igcgp2)=ggitp(:,igcgp2)+ggg(iggg)*geigq1(igp1,:) !time-consuing
      enddo
      enddo
      ggitp = dconjg(ggitp)
!! zmelp <=  \sum_G2 ggitp(Gc+G2) geigqg2(G2))
!! note \bfr'= g (\bfr) +\delta_g  (\bfr= {\bf r})
!! mapping of function g[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
      allocate(gp2phas2(nt0),phase(ngc))
      do igc=1,ngc
        phase(igc)=exp( img*2d0*pi*sum((q+matmul(qlat,ngveccR(:,igc)))*shtv) )
      enddo  
!! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
      allocate(zmelp0(ngc,nt0,ntp0))
      zmelp0=0d0
      allocate(ggitp_(ngp2,ngc)) 
      do itp= 1,ntp0
        do igc=1,ngc
            do igp2=1,ngp2
               nn = ngveccR(:,igc) + ngvecp2(:,igp2)
               igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
               ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
            enddo
         enddo
         zmelp0(:,:,itp)= matmul(transpose(ggitp_),geigq2)
      enddo
      deallocate(ggitp_,ngveccR,phase,ggitp,gp2phas2)
      call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
      deallocate(zmelp0)
      if(verbose()>=100) write(6,*)' melpln2t: end'
      end subroutine melpln2t
!sssssssssssssssssssssssssssssssssssssssss      
      subroutine psi2b_v3(nctot,ncc,nt0,ntp0,iclass,phase,
     i                   cphik, !rbk,cbk,rhbk,chbk,
     i                   cphikq,!rbkq,cbkq,rhbkq,chbkq,
     i                   ppb, !pdb,dpb,ddb,
     i                   mdim,imdim,iatomp,
     d                   mdimx,nlmto,nbloch,nlnmx,
     d                   natom,nclass,
     o zpsi2b)
c originaly 92.03.17 by Ferdi. takao modified at Apr 2002(v2) Feb2006(v3).
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
      integer :: nlmto,nctot,ncc,nt0,ntp0,nlnmx,natom,nbloch,
     &     ia,ic,nc,nv,nc1,ias,iap,icp,i,mdimx,nclass,itp,jp,ib,
     &     nzwork1,nzwork2,mdim(nclass),iclass(natom),
     i     imdim(natom),iatomp(natom)
      integer,allocatable::iasx(:)
      real(8)::  ppb(nlnmx,nlnmx,mdimx,nclass)
      complex(8):: cphik(nlmto,*),cphikq(nlmto,ntp0),phase(natom),
     &     zpsi2b(nbloch,nt0+nctot,ntp0+ncc)
      complex(8),allocatable :: zz(:,:), zwork(:,:),zppb(:,:)
      allocate( zppb(nlnmx,nlnmx),zz(nlnmx,ntp0))
      if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
      if(sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
      allocate(iasx(natom))
      ias = 1
      do ia = 1,natom
        iasx(ia) = ias
        ias = ias + nlnmv(iclass(ia))
      enddo
      if(ias-1/=nlmto) call rx( ' psi2b_v3:sum(nlnmv)/= nlmto')
      do  ia = 1,natom
         ic   = iclass(ia)
         nc   = nlnmc(ic)
         nv   = nlnmv(ic)
         nc1  = nc + 1
         ias  = iasx(ia)
         iap  = iatomp(ia)
         icp  = iclass(iap)
         do   i = 1,mdim(icp)   ! loop over optimal product basis
!     sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)> !bkq is complex but < > is real
           ib = imdim(iap)-1+i
           zpsi2b(ib,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)= phase(ia) * !!     <psi(k+q,t') | psi(k,t) B(i)>
     &        matmul( transpose(cphik(ias:ias+nv-1,1:nt0)), 
     &        dconjg( matmul(transpose(ppb(nc1:nc+nv,nc1:nc+nv,i,icp)),cphikq(ias:ias+nv-1,1:ntp0)) ) ) !zz
         end do                 !end of optimal product basis-loop
      end do !end of atom-loop
      deallocate(zz,zppb,iasx)
      end
c------------------------------------------------------------------------------------
      subroutine psicb_v3 (nctot,ncc,nt0,ntp0,  iclass, phase, 
     i                   cphik,
     i                   cphikq, 
     i                   ppb,
     i                   mdim,
     i                   imdim,iatomp,
     i                   mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i                   icore,ncore,nl,nnc,  
     o            zpsi2b)
      implicit none
      intent(in)          nctot,ncc,nt0,ntp0,  iclass, phase, 
     i                   cphik,
     i                   cphikq, 
     i                   ppb,
     i                   mdim,
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
      integer:: nl,nlmto,nbloch,nclass,natom,nctot,ncc,nt0,ntp0,mdimx,nlnmx,nnc,
     & ib,ias,ics,ia,ic,nc,nv,nc1,iap,icp,i,itp,it,icr,icore(nl*nl*nnc,nclass),ncore(nclass),
     & mdim(nclass),iclass(natom),imdim(natom),iatomp(natom),verbose,ixx
      real(8)::   ppb(nlnmx,nlnmx,mdimx,nclass) !  coskt(natom),sinkt(natom)
      complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc),
     &             cphikq(nlmto,*), cphik(nlmto,*),phase(natom)
      if(ncc/=0 .and. ncc/=nctot) call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
      if(sum(ncore(iclass(1:natom)))/= nctot) call rx( "psicb_v3:sum(ncore) wrong")
      zpsi2b = 0d0
      ias        = 1
      ics        = 0
      do ia = 1,natom
        ic         = iclass(ia)
        nc         = nlnmc(ic)
        nv         = nlnmv(ic)
        nc1        = nc + 1
        iap        = iatomp(ia)
        icp        = iclass(iap)
        ib         = imdim(iap)-1
        do i = 1,mdim(icp)
          ib         = ib + 1
          do itp = 1,ntp0
          do  it = 1,ncore(ic)
              icr      = icore(it,ic) 
              zpsi2b(ib,ics+it,ncc+itp) = phase(ia)* 
     &             dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))
          enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop
          if(ncc==0) cycle
          do itp = 1,ncore(ic)
          do  it = 1,nt0
              icr      = icore(itp,ic)
              zpsi2b(ib,nctot+it,ics+itp) = dconjg(phase(ia))* 
     &             sum(cphik(ias:ias+nv-1,it)*ppb(nc1:nc+nv,icr,i,icp))
          enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop
        enddo !end of optimal product basis-loop
        ias   = ias + nlnmv(ic)
        ics   = ics + ncore(ic)
      enddo
      end
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine drvmelp3( q, ntp0, q_rk,nt0, qik, isp,ginv,
     i        ngc,ngcmx,ngpmx,nband,itp0, 
     i        symope, shtv, qlat, qlatinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,noccx,nctot,
     o       zmelt) 
      use m_readqg,only: readqg
      use m_readeigen,only:readgeigw
! this is for Wanner (readeigW, drvmelp3)      
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
      real(8):: qlat(3,3),shtv(3),qdiff(3),add(3)
     &         ,qibz(3,nqibz),qbz(3,nqbz),qlatinv(3,3),symope(3,3)
     &         ,rmel(nbloch,noccx,ntp0) 
     &         ,cmel(nbloch,noccx,ntp0) 
     & ,pi=3.1415926535897932D0 
      complex(8):: geig1(ngpmx,nband),geig2(ngpmx,nband)
     &            ,zmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      integer(4)::verbose
      call readqg('QGpsi', q,    qt,   ngp1, ngvecpB1)
      call readqg('QGpsi', q_rk, q_rkt,ngp2, ngvecpB2)
      call readgeigW(q,    ngpmx, isp, qu1, geig1)
      call readgeigW(q_rk, ngpmx, isp, qu2, geig2)
      if(sum(abs(qt-qu1))>1d-10) stop 'drvmelp3;qu1/=qu1x'
      if(sum(abs(q_rkt-qu2))>1d-10) stop 'drvmelp3;qu2/=qu2x'
      if(verbose()>=100) write(6,*)' end of read geig '
      qdiff = matmul(symope,qik)  - qt + q_rkt ! rk    -q  +(q-rk) is not zero.
      nadd  = nint(matmul(qlatinv,qdiff))
      allocate( zmelpl(ngc,nt0,ntp0) )
      if (nt0 .ne. ntp0) stop 'drvmelp3: nt0 .ne. ntp0'
      call melpln2t( ngp1, ngvecpB1
     &           ,   ngp2, ngvecpB2 
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,  ! q1=q    ---> iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,itp0(1:nt0)), nt0,     ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
     i       shtv, matmul(symope,qik),qik, symope, qlat, qt, !qt is dummy...
     o       zmelpl)
      zmelt=0d0
      zmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = 
     &  dcmplx(rmel (1:nbloch, 1:nctot+nt0, 1:ntp0),  cmel (1:nbloch, 1:nctot+nt0, 1:ntp0))
      zmelt(nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0) =zmelpl(1:ngc,  1:nt0, 1:ntp0)
      deallocate(zmelpl) 
      end
      
      end module m_zmel
