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
      use m_readQG,only: ngpmx,ngcmx
c      use m_readclasst,only: Readclasst,iclasst !iclass and iclasst is the same?
      use m_hamindex0,only: Readhamindex0,iclasst !iclass and iclasst is the same?
!! q-dependent
      use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
!!------------------------------------------------------
!! SUBROUTINES
      public:: Get_zmel_init, Get_zmel_modex0,  Dconjg_zmel, Deallocate_zmel, 
     &         Setppovlz, Setppovlz_chipm,  Mptauof_zmel ,rwzmel, !GramSchmidt_zmel,
     &   drvmelp3,ppbafp_v2 !for old-fashioned codes in wannier/ 
!! OUTPUT:  zmel for exchange=F, zmeltt for exchange=T.
      complex(8),allocatable,protected,public :: zmel(:,:,:)   !for  Get_zmel
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
      real(8),allocatable,private :: rmelt(:,:,:),cmelt(:,:,:)
!! set by Get_zmel_modex0, which is called right before Get_zmel_init in x0kf_v4h
      logical:: modex0=.false. 
      integer:: nkmin, nkqmin, isp_k, isp_kq,nmtot,nqtot,ispq_bk,ispm_bk
      logical:: debug=.false.

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
      integer:: i!,nbloch
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
      subroutine ppbafp_v2_zmel(ng) !nspin,
      intent(in)::              ng 
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
      subroutine Get_zmel_init(q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc,iprx)
      intent(in)::             q,kvec,irot,rkvec,isp, nmmax,nqmax, nctot,ncc,iprx
!! Get <phiq(q,ncc+nqmax,ispq) |phim(q-rkvec,nctot+nmmax,ispm) MPB(rkvec,ngb)> ZO^-1
!! Right after call get_zmel_init, we get zmel defined in this module.
!! kvec is in the IBZ, rk = Rot_irot(kvec)
!! \parameter all inputs
!! \parameter matrix <MPB psi|psi>
      logical:: iprx
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
      if(iprx) then
         print *,'ppp modex0',modex0
         write(6,"(a,l,10i4)")'pppp3:',irot,nmini,nmmax,ispm,nctot,nqini,nqmax,ispq,ncc
         write(6,"(a,10(3f9.3,x))") 'pppp3:',q, kvec,rkvec
      endif         
      call get_zmelt2(
     &    kvec=kvec,irot=irot,rkvec=rkvec,          ! MPB_rkvec
     &    nmini=nmini,nmmax=nmmax,ispm=ispm ,nctot=nctot,  ! middle-phi for  phi_{q-rkvec}
     &    q=q,nqini=nqini,nqmax=nqmax,ispq=ispq,ncc=ncc  )! end-phi    for  phi_q
      end subroutine get_zmel_init
      
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
c      logical:: exchange
!! memo =========================================
!! zmelt = rmelt(igb(rkvec), iocc(q), iunocc(q-rkvec)) + i* cmelt
!! iunocc: band index at target  q.
!! iocc:   band index at intermediate vector qk = q - rkvec 
!! igb: index of mixed product basis       at rkvec (or written as rk)
!!   igb=1,ngb 
!!   ngb=nbloch+ngc  ngb: # of mixed product basis
!!                   nbloch: # of product basis (within MTs)
!!                   ngc: # of IPW for the Screened Coulomb interaction.
!!                   igc is for given 
! zmeltt= <itp|it,ib>=(rmelt,-cmelt)
! zmel(igb,it,itp) = C(ppovlz)*N(zmeltt(:,it*itp)) ! C means Hermitian conjugate, N means normal
      subroutine get_zmelt2(
     &   kvec,irot,rkvec,          ! rkvec MPB  ngcxxx,ngbxxx,  MPB_rkvec
     &   nmini,nmmax,ispm,nctot,   ! q-rkvec middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq,ncc)  ! end state for phi_q
      intent(in)::          
     &   kvec,irot,rkvec,          
     &   nmini,nmmax,ispm,nctot,   
     &    q,nqini,nqmax,ispq,ncc   
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
      if(debug) write(*,*) 'get_zmelt2 in m_zmel: start'
      if(.not.done_genallcf_v3) call rx('m_zmel: not yet call genallcf_v3')
      if(.not.done_rdpp)        call rx('m_zmel: not yet call rdpp')
      if(.not.done_read_bzdata) call rx('m_zmel: not yet call read_bzdata')

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
c      if(debug) write(6,'("Goto psi2b_v3 nctot ncc nmmax nqmax=",4i4)') nctot,ncc,nmmax,nqmax
c      if(debug) write(6,'("4444 zzzmelsum ",3i5,3d13.5)') nbloch,nmtot,nqtot,sum(abs(zzzmel)),sum(zzzmel)
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
c      if(debug) write(6,'("5555 zzzmelsum ",3i5,3d13.5)')nbloch,nmtot,nqtot,sum(abs(zzzmel)),sum(zzzmel)
c      if(debug) write(6,'("6666 chim      ",4i5,4d13.5)')nmini,nmmax,nqini,nqmax,sum(cphim(:,nmini:nmmax)),sum(cphiq(:,nqini:nqmax))
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
     i  symope, shtv, qlat, qlatinv,qibz,qbz,nqbz,nqibz,
     i  drealzzzmel, dimagzzzmel, nbloch, nctot,ncc,itps,
     o  rmelt,cmelt)
      deallocate(drealzzzmel,dimagzzzmel)
      if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)
      allocate( zmel (ngb, nmtot, nqtot) )
      call matmc(ppovlz, dcmplx(rmelt,-cmelt), zmel,ngb,ngb,nmtot*nqtot)
      deallocate(rmelt,cmelt)
      end subroutine

      
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
c      integer(4):: zvztest
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
c                  if(zvztest()==2) then
c                    if(abs(ppb(i1,i2,i,ic))/=0d0  )then
c                      write(6,"('i1 i2 i ic=',4i4,' ppb=',f10.6)") i1,i2,i,ic, ppb(i1,i2,i,ic)
c                    endif
c                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      end
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
!> Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
!! q1=q    ---> iq 1:ntp0 q-point
!! q2=q-rk ---> kp 1:nt0  occupied
!!  zmelpl
      subroutine drvmelp( q, ntp0, q_rk,nt0, qik, ispq,ispq_rk,ginv,
     i        ngc,ngcmx,ngpmx_in,nband,itp0, 
     i        symope, shtv, qlat, qlatinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,nctot,ncc,itps,
     o       rmelt,cmelt)
      use m_readqg, only : readqg
      use m_readeigen,only : readgeigf
c-----------------------------------------------------
      implicit none
      real(8),intent(in):: q(3),q_rk(3),qik(3),ginv(3,3)
      integer(4),intent(in):: ngpmx_in,nqbz,nqibz,ngcmx ,nctot,nband,
     &  ntp0,nt0,nbloch,ngc,itp0(ntp0),ispq,ispq_rk,itps,ncc
      real(8),intent(in):: qlat(3,3),shtv(3),qibz(3,nqibz),qbz(3,nqbz),
     &  qlatinv(3,3),symope(3,3),rmel(nbloch,nctot+nt0,ncc+ntp0),
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
      geig1= readgeigf(q,ispq)       !call readgeig(q,    ngpmx_in, ispq, qu1, geig1)
      geig2= readgeigf(q_rk,ispq_rk) !call readgeig(q_rk, ngpmx_in, ispq_rk, qu2, geig2)
      qdiff = matmul(symope,qik)  - qt + q_rkt
! rk    -q  +(q-rk) is not zero. <rk q-rk |q>
      nadd = nint(matmul(qlatinv,qdiff))
c      if(sum(abs(add-nadd))>tolq) call rx('drvmelp add/=nadd')
      if(ngc/=0) then
        allocate( zmelpl(ngc,nt0,ntp0) )
        call melpln2t(ngp1, ngvecpB1
     &   ,  ngp2, ngvecpB2
     &   ,   ngc,  nadd,
     &   geig1(1:ngp1,itps-1+itp0(1:ntp0)), ntp0, ! q1=(shifted q) ->iq ngp1 1:ntp0 q-point
     &   geig2(1:ngp2,1:nt0), nt0, ! q2=(shifted q-rk) -> kp ngp2 1:nt0  occupied 
     i   shtv, matmul(symope,qik),qik, symope, qlat,
     i   qt,                    !qt oct2013 means q1
     o   zmelpl)
      endif
      rmelt=0d0
      cmelt=0d0
      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
c      if(ngc/=0) then
        rmelt(nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)= dreal( zmelpl(1:ngc, 1:nt0,  1:ntp0) )
        cmelt(nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)= dimag( zmelpl(1:ngc, 1:nt0,  1:ntp0) )
        deallocate(zmelpl)
c      endif
c      if(debug) write(6,*)'end of drvmelp'
      end

!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
!> Mattrix elements <Plane psi |psi> from interstitial plane wave.
!! zmelp(igc(qi),it(q2),itp(q1)) = <itp(for q1+G1)| it(for q2+G2) igc> 
!! NOTE: shtv = g(delta_{g^-1})
      subroutine melpln2t
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied                        
     i       shtv,q, qi, symope, qlat, !
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
      real(8),intent(in) :: shtv(3),q(3),qi(3), symope(3,3),qlat(3,3)
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
!!      call getppx(qlat,  ngcx,qi,symope) ! rppovl.f

!! allocate and return foobar decleared in use m_read_ppovl,only: foobar
      if(debug) write(6,*)' melpln2t: nadd=',nadd
      call getppx2(qlat,qi) ! rppovl.f
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
      call rotgvec(symope, 1, ngc, ngcs, qlat, ngvecc, 
     o    ngveccR) 
      if(verbose()>=100) write(6,*)' melpln2t: godo igp1 loop'
      allocate(ggitp(ntp0,ngcgp2))
      ggitp = 0d0
      do igp1  = 1,ngp1   !for ngp1
      do igcgp2= 1,ngcgp2 !for ngc+ngp2 
         nn = ngvecp1(:,igp1)- nvgcgp2(:,igcgp2) - nadd 
         ! G1 -(Gc+G2) - Gadd !  -Gadd= -rk + qt -q_rk
         if(nn(1)<nxminggg .or. nxmaxggg<nn(1) ) cycle
         if(nn(2)<nyminggg .or. nymaxggg<nn(2) ) cycle
         if(nn(3)<nzminggg .or. nzmaxggg<nn(3) ) cycle
         iggg = igggi(nn(1),nn(2),nn(3))
         if(iggg<0) cycle
!! ggg(iggg) = <qt+G1 -(rk+Gc) -(q_rk+G2) >, where 
!!       gggx(1:ntp0,igcg2) = gggx(1:ntp0,igcg2)+ ggg(iggg)*geigq1_(1:ntp0,igp1)
! this is time-consuging part
         call zaxpy(ntp0,ggg(iggg), geigq1(igp1,1),ngp1, ggitp(1,igcgp2),1)
      enddo
      enddo
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
        phase(igc)=exp( img*2d0*pi*sum((q+matmul(qlat,ngveccR(:,igc)))*shtv) )
      enddo  
!! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
      allocate(zmelp0(ngc,nt0,ntp0))
      zmelp0=0d0
!! If range of nvgcgp2 is not large enough, igcgp2 can be -100000.
      allocate(z2(ngc,nt0))
      allocate(ggitp_(ngp2,ngc)) !gp2))ngcgp2))
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
      deallocate(ggitp_)
      if(verbose()>=100) write(6,*)' melpln2t: 4444444444'
      deallocate(ngveccR,phase)
      call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
      deallocate(igcgp2i,ggitp,gp2phas2,zmelp0,z2)
      if(verbose()>=100) write(6,*)' melpln2t: end'
      end subroutine melpln2t
!sssssssssssssssssssssssssssssssssssssssss      
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
      allocate( zppb(nlnmx,nlnmx),zz(nlnmx,ntp0))
      if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
      if( sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
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
         if (nc+ nlnmv(ic) > nlnmx) call rx( 'psi2b_v3: nlnmx exceeded')
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
        call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
      endif
      if(sum(ncore(iclass(1:natom)))  /= nctot) 
     &   call rx( "psicb_v3:sum(ncore) wrong")
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
        iap        = iatomp(ia)
        icp        = iclass(iap)
        ib         = imdim(iap)-1
        do       i = 1,mdim(icp)
          ib         = ib + 1
          do     itp = 1,ntp0
          do      it = 1,ncore(ic)
              icr      = icore(it,ic) ! write(6,*),it,ic,icore(it,ic)
              zpsi2b(ib,ics+it,ncc+itp) = phase(ia)* 
     &             dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))
          enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop
          if(ncc==0) cycle
          do     itp = 1,ncore(ic)
            do      it = 1,nt0
              icr      = icore(itp,ic) ! write(6,*),it,ic,icore(it,ic)
              if(verbose()>90) write(6,*)"psicb_v3: aaa itp it ic icr=",itp,it,ic,icr
              zpsi2b(ib,nctot+it,ics+itp) = dconjg(phase(ia))* 
     &             sum(cphik(ias:ias+nv-1,it)*ppb(nc1:nc+nv,icr,i,icp))
            enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop
          if(verbose()>90) write(6,*)"psicb_v3: xxx ia i=",ia,i
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
     & dcmplx(rmel (1:nbloch, 1:nctot+nt0, 1:ntp0),  cmel (1:nbloch, 1:nctot+nt0, 1:ntp0))
      zmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0) =zmelpl(1:ngc,  1:nt0, 1:ntp0)
      deallocate(zmelpl) 
      end
      
      end module m_zmel
