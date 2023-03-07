!! Get the matrix element zmel =  ZO^-1 <MPB psi|psi> , where ZO is ppovlz.
!!  "call get_zmel_init" return zmel
!!  All dependencies (use foobar below ) are inputs (must be protected).
module m_zmel
  use m_genallcf_v3,only: &
       nclass,natom,nspin,nl,nn,nnv,nnc, &
       nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw, &
       alat,delta,deltaw,esmr,iclass,nlnmv,  nlnmc,  icore,ncore,occv,unoccv , &
       occc,unoccc, nocc, nunocc, plat, pos,z,ecore, &
       il, in, im, mnl=>nlnm , nl,nn,nlnmx
  use m_hamindex,only: ngrp, symgg=>symops,invg=>invgx
  use m_rdpp,only: Rdpp, nxx,lx,nx,mdimx,nbloch,cgr,ppbrd,nblocha,done_rdpp
  use m_readeigen,only: Readcphif 
  use m_read_bzdata,only: nqbz,nqibz,  qlat,ginv,qbz,qibz,wbz, done_read_bzdata
  use m_readhbe,only: nband
  use m_itq,only: itq,ntq
  use m_readQG,only: ngpmx,ngcmx,Readqg
  use m_hamindex0,only: Readhamindex0,iclasst
  use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
  ! OUTPUT: zmel(nbb,nmtot, nqtot) ,nbb:mixproductbasis, nmtot:middlestate, nqtot:endstate
  complex(8),allocatable,protected,public :: zmel(:,:,:),zmel0(:,:,:) !for  Get_zmel
  !                                                      zmel0 is just one another zmel which can be contained in m_zmel
  real(8),protected,public:: qm0(3) !for zmel0
  integer,protected:: nbb           !1st dimension of zmel. MPB
  public:: Get_zmel_init, Dconjg_zmel, Deallocate_zmel, Deallocate_zmel0, &
       Setppovlz, Setppovlz_chipm,  Mptauof_zmel,rwzmel &! & GramSchmidt_zmel (test code)
       ,&
       drvmelp3,setzmel0,unsetzmel0,ppbafp_v2 !these are for old-fashioned codes in wannier/
  !!------------------------------------------------------
  !! set by Mptauof_zmel in advance
  private
  integer,allocatable,private :: miat(:,:)
  real(8),allocatable,private :: tiat(:,:,:),shtvg(:,:)
  real(8),allocatable,private :: ppbir(:,:,:,:,:,:)
  complex(8),allocatable,private :: ppovlz(:,:)
  real(8),private:: qlatinv(3,3),q_bk(3)=1d10,qk_bk(3)=1d10
  logical,private:: init=.true.
  complex(8),allocatable,private :: cphiq(:,:), cphim(:,:),cphitemp(:,:)
!  logical:: modex0=.false.
  integer:: nkmin, nkqmin, isp_k, isp_kq,nmtot,nqtot,ispq_bk,ispm_bk
  logical:: debug=.false.,zzmel0=.false.
  integer:: nbbx=0
contains
  ! ssssssssssssssssssssssssssssssssss
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
  end subroutine GramSchmidt_zmel
  subroutine Dconjg_zmel()
    if(zzmel0) then
       zmel0 = dconjg(zmel0)
    else
       zmel = dconjg(zmel)
    endif
  end subroutine Dconjg_zmel
  subroutine Deallocate_zmel()
    if(allocated(zmel)) deallocate(zmel)
  end subroutine Deallocate_zmel
  subroutine Deallocate_zmel0()
    if(allocated(zmel0)) deallocate(zmel0)
  end subroutine Deallocate_zmel0
  subroutine setppovlz(q,matz)
    intent(in)::       q,matz
    logical:: matz
    !! set ppovlz for given q
    !    ppolvz(igb,ivcou)= (1    0 ) \times  zcousq(igb, ivcou)
    !                       (0 ppovl)
    !    If matz=F, no multiplication by ivcou.  Thus we have ppolz(igb,igb)
    real(8) :: q(3)
    complex(8),allocatable :: ppovl_(:,:),ppovl(:,:)!,ppovlzinv(:,:)
!    logical:: eibz4x0
    integer:: i
    if(allocated(ppovlz)) deallocate(ppovlz)
    if(allocated(ppovl)) deallocate(ppovl)
    allocate( ppovl(ngc,ngc),ppovlz(ngb,ngb))
    call readppovl0(q,ngc,ppovl) !q was qq
    ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
    ppovlz(nbloch+1:nbloch+ngc,:)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
    deallocate(ppovl)
  end subroutine setppovlz
  subroutine setppovlz_chipm(zzr,nmbas1)
    intent(in)::             zzr,nmbas1
    integer::nmbas1
    complex(8):: zzr(ngb,nmbas1)
    if(allocated(ppovlz)) deallocate(ppovlz)
    allocate(ppovlz(ngb,nmbas1))
    ppovlz= zzr
    nbbx=nmbas1
  end subroutine setppovlz_chipm
  subroutine mptauof_zmel(symops,ng)! Set miat,tiat,invgx,shtvg, and then call ppbafp_v2_zmel successively
    intent(in)::          symops,ng
    integer:: ng
    real(8):: symops(9,ng)
    integer,allocatable ::  invgx(:)
    call readhamindex0()
    allocate(invgx(ng),miat(natom,ng),tiat(3,natom,ng),shtvg(3,ng))
    call mptauof(symops,ng,plat,natom,pos,iclasst &
         ,miat,tiat,invgx,shtvg ) !note: miat,tiat,shtvg are defined in m_zmel.
    deallocate(invgx)
    call Rdpp(ng,symops)
    call ppbafp_v2_zmel(ng)
  end subroutine mptauof_zmel
  !----------------------------------------------------
  subroutine ppbafp_v2_zmel(ng)
    intent(in)::            ng
    integer :: ng,is,irot,lmxax
    integer :: ic, i,lb,nb,mb,lmb,i1,ibas,i2 !nl,nn,
    integer :: np,lp,mp,lmp,n,l,m,lm
    lmxax=nl-1
    allocate( ppbir(nlnmx,nlnmx,mdimx,nclass,ng,nspin))
    do irot = 1,ng
       do is = 1,nspin 
          !! ppbir is rotated <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)> by rotated cg coefficients cgr
          do concurrent (ic=1:nclass)
             ibas = ic
             i = 0 !i = product basis index.
             do lb  = 0, lx (ibas)
                do nb  = 1, nx (lb,ibas)
                   do mb  = -lb, lb
                      i    = i+1           !The number of product basis is  =(i at the end of loop).
                      lmb  = lb*lb + lb + mb + 1
                      do concurrent (i2 = 1:mnl(ic),i1 = 1:mnl(ic)) !phi1 phi2 index
                         np   = in(i2,ic)
                         lp   = il(i2,ic)
                         mp   = im(i2,ic)
                         lmp  = lp*lp + lp + mp + 1
                         n    = in(i1,ic)
                         l    = il(i1,ic)
                         m    = im(i1,ic)
                         lm   = l*l + l + m + 1
                         ppbir(i1,i2,i,ic,irot,is)=cgr(lm,lmp,lmb,irot)*ppbrd(l,n,lp,np,lb,nb,is+nspin*(ic-1))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ppbafp_v2_zmel
  subroutine ppbafp_v2 (irot,ng,isp, mdimx,lx,nx,nxx, &! & Bloch wave
    cgr,lmxax,    & !rotated CG
    ppbrd,        & !radial integrals
    ppb)
    ! <Phi(RLn) Phi(RL'n') B(R,i)>. Atomic part of MPB
    !  n differenciates core phi phidot localOrbital.
    !  in,il,im      = index for n,l,m s. indxlnm.f
    !   ppb            = <Phi(RLn) Phi(RL'n') B(R,i)>
    implicit none
    integer,intent(in) :: irot,ng,isp!,nspin,nclass,nlnmx,mdimx
    integer,intent(in) :: lx(nclass),nx(0: 2*(nl-1),nclass)
    integer,intent(in) :: lmxax,nxx,mdimx
    real(8), intent(in) :: cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng)
    real(8), intent(in) :: ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nxx,nclass*nspin)
    real(8), intent(out) :: ppb(nlnmx,nlnmx,mdimx,nclass)
    integer :: ic, i,lb,nb,mb,lmb,i1,ibas,i2, np,lp,mp,lmp,n,l,m,lm
    do concurrent (ic=1:nclass)
       ibas = ic
       i = 0 !i = product basis index.
       do lb  = 0, lx (ibas)
          do nb  = 1, nx (lb,ibas)
             do mb  = -lb, lb
                i    = i+1           !The number of product basis is  =(i at the end of loop).
                lmb  = lb*lb + lb + mb + 1
                do concurrent (i2 = 1:mnl(ic),i1 = 1:mnl(ic)) !phi1 phi2 index
                   np   = in(i2,ic)
                   lp   = il(i2,ic)
                   mp   = im(i2,ic)
                   lmp  = lp*lp + lp + mp + 1
                   n    = in(i1,ic)
                   l    = il(i1,ic)
                   m    = im(i1,ic)
                   lm   = l*l + l + m + 1
                   ppb(i1,i2,i,ic)=cgr(lm,lmp,lmb,irot)*ppbrd(l,n,lp,np,lb,nb,isp+nspin*(ic-1))
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ppbafp_v2
  subroutine rwzmel(iq,k,isp_k,isp_kq,rw,q,izmel0)
    intent(in)::    iq,k,isp_k,isp_kq,rw,q,izmel0
    logical,optional::izmel0
    integer:: iq,k,isp_k,isp_kq,izmel,nmtot_,nqtot_
    character(1):: rw
    character(10) :: i2char
    character*256::fff,fname
    real(8),optional:: q(3)
    integer::nnn
    fname=trim(i2char(iq))//'_'//trim(i2char(k))//'_'//trim(i2char(isp_k))//trim(i2char(isp_kq))
    open(newunit=izmel,file='zmel.'//trim(fname),form='unformatted')
    if(rw=='r' .AND. present(izmel0)) then
       read(izmel) nbb,nmtot_,nqtot_,qm0
       if(allocated(zmel0)) deallocate(zmel0)
       allocate(zmel0(1:nbb,1:nmtot_,1:nqtot_))
       read(izmel) zmel0
    elseif(rw=='r') then
       read(izmel) nbb,nmtot_,nqtot_,qm0
       if(allocated(zmel)) deallocate(zmel)
       allocate(zmel(1:nbb,1:nmtot_,1:nqtot_))
       read(izmel) zmel
    else
       write(izmel) nbb,nmtot,nqtot,q
       write(izmel) zmel
    endif
    close(izmel)
  end subroutine rwzmel
  subroutine Get_zmel_init(q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc,iprx)
    intent(in)::           q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc,iprx
    ! ns1:ns2 is the range of middle states (count including core index (i=1,nctot). Thus kth valence is at nctot+k.
    !! Get zmel= <phiq(q,ncc+nqmax,ispq) |phim(q-rkvec, ns1:ns2, ispm) MPB(rkvec,ngb)> ZO^-1
    !! kvec is in the IBZ, rk = Rot_irot(kvec)
    !! \parameter all inputs
    !! \parameter matrix <MPB psi|psi>
    !! igb: index of mixed product basis       at rkvec (or written as rk)
    !!   igb=1,ngb
    !!   ngb=nbloch+ngc  ngb: # of mixed product basis
    !!                   nbloch: # of product basis (within MTs)
    !!                   ngc: # of IPW for the Screened Coulomb interaction.
    !!                   igc is for given
    ! zmelt= <itp|it,ib>
    ! zmel(igb,it,itp) = transpose(dcojng(ppovlz))*zmelt(:,it,itp)
    !   matmul(symgg(:,:,irot),kvec)-rkvec can have difference of reciprocal vectors.
    logical:: iprx
    integer:: isp,nmmax,ns1,ns2,nqmax,irot,ispq,ispm,nmini,nqini, nctot,ncc
    real(8) ::  quu(3),q(3), kvec(3),rkvec(3)
    real(8),parameter::tolq=1d-8
    nmini = ns1          !starting index of middle state  (nctot+nvalence order)
    nmmax = ns2-nctot    !end      index of middle state  (nctot+nvalence order)
    ZmeltMain: Block
      use m_readeigen,only : readgeigf
      integer:: it,ia,kx,verbose,nstate,imdim(natom),nt0,ntp0,invr, iatomp(natom),ispq_rk
      real(8):: symope(3,3),shtv(3),tr(3,natom),qk(3)
      real(8),allocatable :: ppb(:,:,:,:)
      complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
      complex(8):: expikt(natom)
      complex(8),allocatable::  zmelt(:,:,:)
      integer :: ngp1, ngp2
      integer :: ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
      real(8):: q_rkt(3),qt(3),q_rk(3)!,qik(3)
      real(8) :: qdiff(3)
      complex(8) :: geig1(ngpmx,nband),geig2(ngpmx,nband)
      logical:: debug=.false.
      if(init) then
         call minv33(qlat,qlatinv)
         allocate( cphiq(nlmto,nband), cphim(nlmto,nband), cphitemp(nlmto,nband))
         init=.false.
      endif
      if(sum(abs(q-q_bk))>tolq .OR. ispq/=ispq_bk)  then
         cphitemp= readcphif(q,ispq)
         cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq))
         q_bk=q
         ispq_bk=ispq
      endif
      qk =  q - rkvec ! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      if(sum(abs(qk-qk_bk))>tolq .OR. ispm/=ispm_bk) then
         cphim= readcphif(qk, ispm) !all readin but we need only bands nmini:nmmax
         qk_bk= qk
         ispm_bk= ispm
      endif
      invr  = invg(irot)       !invrot (irot,invg,ngrp) ! Rotate atomic positions invrot*R = R' + T
      tr    = tiat(:,:,invr)
      iatomp= miat(:,invr)
      symope= symgg(:,:,irot)
      shtv  = matmul(symope,shtvg(:,invr))
      allocate( ppb(nlnmx,nlnmx,mdimx,nclass)) ! ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
      ppb = ppbir(:,:,:,:,irot,ispq) !rotated ppb. MPB has no spin dependence
      do ia = 1,natom
         imdim(ia)  = sum(nblocha(iclass(1:ia-1)))+1
         expikt(ia) = exp(img *tpi* sum(kvec*tr(:,ia)) )  !phase expikt(ia) is for exp(ik.T(R))
      end do
      nt0= nmmax-nmini+1
      ntp0=nqmax-nqini+1
      nmtot  = nctot + nt0     ! = phi_middle nmtot=ns2-ns1+1
      nqtot  = ncc   + ntp0    ! = phi_end
      allocate(zmelt(1:nbloch+ngc,nmtot,nqtot))
      zmelt=0d0
      CoreRelatedPartWithinMT: if(ncc>0 .OR. nctot>0) then
         call psicb_v3  ( nctot,ncc,nmmax,nqmax,iclass,expikt, &
              cphim(1,nmini), & ! middle phi
              cphiq(1,nqini), & ! end phi
              ppb,  nblocha, imdim,iatomp, &
              mdimx,nlmto,nbloch,nlnmx,natom,nclass, &
              icore,ncore,nl,nnc, &
              zmelt(1:nbloch,:,:))
      endif CoreRelatedPartWithinMT
      ValenceValenceWithinMT: if(nt0*ntp0>0) then ! nt0=valence num of middle states ! ntp0=valence num of end states
         call psi2b_v3( nctot,ncc, nt0, ntp0, iclass,expikt, &
              cphim(1,nmini), cphiq(1,nqini), ppb, nblocha, imdim,iatomp, mdimx,nlmto,nbloch,nlnmx, natom,nclass, &
              zmelt(1:nbloch,:,:))
         ! IPW part for Valence x Valence = <IPW phi_valence |phi_valence>
         q_rk=qk
         ispq_rk=ispm
         call readqg('QGpsi',q,    qt,   ngp1, ngvecpB1) !q is mapped to qt in BZ
         call readqg('QGpsi',q_rk, q_rkt,ngp2, ngvecpB2)
         geig1= readgeigf(q,ispq)       !call readgeig(q,    ngpmx_in, ispq, qu1, geig1)
         geig2= readgeigf(q_rk,ispq_rk)
         qdiff = matmul(symope,kvec)  - qt + q_rkt ! rk -q +(q-rk) is not zero. <rk q-rk |q>
         nadd = nint(matmul(qlatinv,qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
         call melpln2t(ngp1, ngvecpB1 &
              ,  ngp2, ngvecpB2,   ngc,  nadd, &
                                ! bugfix (probably affects to interband,intraband dielectric function cases?) 2022-8-28
                                !           geig1(1:ngp1,nqini-1+itq(1:ntp0)), ntp0,& ! &  q1=(shifted q) ->iq ngp1 1:ntp0 q-point
                                !           geig2(1:ngp2,1:nt0), nt0, &! &  q2=(shifted q-rk) -> kp ngp2 1:nt0  occupied
              geig1(1:ngp1, itq(nqini:nqmax)),ntp0,& ! &  q1=(shifted q) ->iq ngp1 1:ntp0 q-point
              geig2(1:ngp2,     nmini:nmmax), nt0, & ! &  q2=(shifted q-rk) -> kp ngp2 1:nt0  occupied
              shtv, matmul(symope,kvec),kvec, symope, qlat, qt, &
              zmelt(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0))
      endif ValenceValenceWithinMT
      nbb=ngb
      if(nbbx/=0) nbb=nbbx
      if(allocated(zmel)) deallocate(zmel)
      allocate(zmel(nbb,ns1:ns2, nqtot))
      if(nt0==0) then ! core only for middle states. fast mode 
         zmel=0d0
         block
           integer:: ia,it,ics,iap,icp
           do concurrent(ia=1:natom)
              ics = sum(ncore(iclass(1:ia-1)))
              iap = iatomp(ia)
              icp = iclass(iap)
              do concurrent(it=1:ncore(iclass(ia)))
                 if(ics+it<ns1.or.ns2<ics+it) cycle
                 associate(nb1 => imdim(iap),nb2 => imdim(iap)+nblocha(icp)-1)
                   zmel(1:nbb,ics+it,1:nqtot) = dconjg(matmul(transpose(ppovlz(nb1:nb2,1:nbb)),zmelt(nb1:nb2,ics+it,1:nqtot)))
                 endassociate
              enddo
           enddo
         endblock
      else
         call matm(dconjg(transpose(ppovlz)), dconjg(zmelt), zmel,nbb,ngb,nmtot*nqtot)
      endif
      if(zzmel0) call move_alloc(from=zmel,to=zmel0)
      deallocate(zmelt)
    EndBlock ZmeltMain
  end subroutine get_zmel_init
  subroutine setzmel0()
    zzmel0=.true.
  end subroutine setzmel0
  subroutine unsetzmel0()
    zzmel0=.false.
  end subroutine unsetzmel0
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
  subroutine melpln2t( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd, &
       geigq1, ntp0, &!  q1=q    ---> iq 1:ntp0 q-point
       geigq2, nt0,  &!  q2=q-rk ---> kp 1:nt0  occupied
       shtv,q, qi, symope, qlat, qt, &
       zmelp)
    use m_read_ppovl,only: getppx2,igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze, &
         nvggg,nvgcgp2,ngvecc, nggg,ngcgp,ngcread, ggg,ppovlinv, ngc2,ngvecc2
    implicit none
    integer:: itp,igc
    integer,intent(in) :: ngp1, ngvecp1(3,ngp1), ngp2
    integer,intent(in) :: ngvecp2(3,ngp2), ngc,nadd(3),ntp0,nt0
    complex(8),intent(in) :: geigq1(ngp1,ntp0),geigq2(ngp2,nt0)
    real(8),intent(in) :: shtv(3),q(3),qi(3), symope(3,3),qlat(3,3)
    real(8),intent(in) :: qt(3)
    complex(8),intent(out) :: zmelp(ngc,nt0,ntp0)
    complex(8),parameter :: img=(0d0,1d0)
    real(8), parameter :: pi=3.1415926535897932D0
    integer :: nn(1:3)
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
    if(verbose()>=90) debug= .TRUE. 
    call getppx2(qlat,qi) ! read and allocate PPOVL*
    if(ngc/=ngcread)call rx( 'melpln2: ngc/= ngcx by getppx:PPOVLG')!  write(6,*)qi,ngcread,ngc
    allocate(ngveccR(1:3,1:ngc))
    ngcs(1) = ngc
    call rotgvec(symope, 1, ngc, ngcs, qlat, ngvecc, &
         ngveccR)
    allocate(ggitp(ntp0,ngcgp))
    ggitp = 0d0
    do concurrent (igcgp2=1:ngcgp) !for ngc+ngp2
       do igp1  = 1,ngp1   !for ngp1
          nn = ngvecp1(:,igp1)- nvgcgp2(:,igcgp2) - nadd ! G1 -(Gc+G2) - Gadd ! -Gadd= -rk + qt -q_rk
          if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
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
    phase(:)=[(exp( img*2d0*pi*sum((q+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]
    !! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
    allocate(zmelp0(ngc,nt0,ntp0))
    zmelp0=0d0
    allocate(ggitp_(ngp2,ngc))
    do concurrent (itp= 1:ntp0)
       do concurrent (igc=1:ngc,igp2=1:ngp2)
          nn = ngveccR(:,igc) + ngvecp2(:,igp2)
          igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
          ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
       enddo
       zmelp0(:,:,itp)= matmul(transpose(ggitp_),geigq2)
    enddo
    deallocate(ggitp_,ngveccR,phase,ggitp,gp2phas2)
    call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
    deallocate(zmelp0)
    if(verbose()>=100) write(6,*)' melpln2t: end'
  end subroutine melpln2t
  subroutine psi2b_v3(nctot,ncc,nt0,ntp0,iclass,phase, &
       cphik,cphikq, ppb, mdim,imdim,iatomp, mdimx,nlmto,nbloch,nlnmx, natom,nclass,&
       zpsi2b)
    use m_lgunit,only:stdo,stdl
    ! originaly 92.03.17 by Ferdi. takao modified at Apr 2002(v2) Feb2006(v3).
    ! calculates <psi(k',t') | psi(k,t) B(R,i)>
    ! for all R
    ! psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
    ! B(R,i)   = Bloch orthonormal product basis for atom R
    ! psi(k,t) is stored after nctot
    ! nt0        = no. of valence middle states
    ! ntp0       = no. of valence end states
    ! cphik  b(k) ! middle states
    ! cphikq b(k')!end states
    ! ppb        = <phi(RLn) phi(RL'n') B(R,i)>
    ! nlnmv      = number of l,n,m for valence
    ! nlnmc      = number of n,l,m for core states
    ! mdim       = number of optimal product basis functions for each class
    ! nctot      = total no. allowed core states
    ! nbloch     = total no. optimal product basis
    ! nlnmx      = maximum number of l,n,m
    ! zpsi2b     =  the matrix elements
    !----------------------------------------------------------
    implicit none
    integer :: nlmto,nctot,ncc,nt0,ntp0,nlnmx,natom,nbloch, &
         ia,ic,nc,nv,nc1,ias,iap,icp,i,mdimx,nclass,itp,jp,ib, &
         nzwork1,nzwork2,mdim(nclass),iclass(natom), &
         imdim(natom),iatomp(natom),ierr
    integer,allocatable::iasx(:)
    real(8)::  ppb(nlnmx,nlnmx,mdimx,nclass)
    complex(8):: cphik(nlmto,*),cphikq(nlmto,ntp0),phase(natom), &
         zpsi2b(nbloch,nt0+nctot,ntp0+ncc)
    complex(8),allocatable :: zz(:,:), zwork(:,:),zppb(:,:)
    allocate( zppb(nlnmx,nlnmx),zz(nlnmx,ntp0))
    if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
    if(sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
    allocate(iasx(natom))
    if(sum(nlnmv(iclass(1:natom)))/=nlmto) call rx( ' psi2b_v3:sum(nlnmv)/= nlmto')
    iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)]
    iatomloop: do concurrent(ia = 1:natom)
       ic   = iclass(ia)
       nc   = nlnmc(ic)
       nv   = nlnmv(ic)
       nc1  = nc + 1
       ias  = iasx(ia)
       iap  = iatomp(ia)
       icp  = iclass(iap)
       prodloop: do concurrent(i=1:mdim(icp))   ! loop over optimal product basis
          !     sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)> !bkq is complex but < > is real
          ib = imdim(iap)-1+i  !   <psi(k+q,t') | psi(k,t) B(i)>
          zpsi2b(ib,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)= phase(ia) * matmul( transpose(cphik(ias:ias+nv-1,1:nt0)), &
               dconjg( matmul(transpose(ppb(nc1:nc+nv,nc1:nc+nv,i,icp)),cphikq(ias:ias+nv-1,1:ntp0)) ) ) !zz
       enddo prodloop 
    enddo iatomloop
    deallocate(zz,zppb)
    call flush(stdo)
  end subroutine psi2b_v3
  !------------------------------------------------------------------------------------
  subroutine psicb_v3(nctot,ncc,nt0,ntp0,iclass,phase,&
       cphik,cphikq,ppb,mdim,imdim,iatomp,mdimx,nlmto,nbloch,nlnmx,natom,nclass,icore,ncore,nl,nnc, &
       zpsi2b)
    implicit none
    intent(in)::       nctot,ncc,nt0,ntp0,iclass,phase,&
         cphik,cphikq,ppb,mdim,imdim,iatomp,mdimx,nlmto,nbloch,nlnmx,natom,nclass,icore,ncore,nl,nnc
    intent(out):: zpsi2b
    !- Calculates <psi (k+q,t) |core(k,t) B(R,ibloch)> ---------------------
    !      also   <core(k+q,t) |psi (k,t) B(R,ibloch)>
    !r B(R,i)   = Mixed basis.
    !r core(k,t)= core states
    !r <psi(k+q,t') | core(k,t) B(R,i)> = S[RLn]  cphik(RLn,k+q,t')^*  * ppb
    !r
    !i nt0        = no. k states
    !i ntp0       = no. kq states
    !i cphik coeefficients of MT part or valence eigenfunction.
    !i icore      = index for core states
    !i ncore      = no. core states in each class
    !i ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>
    !i nlnmv      = number of l,n,m for valence
    !i nlnmc      = number of l,n,m for core states
    !i mdim       = number of optimal product basis functions for each class
    !i nbloch     = total no. optimal product basis
    !i nlnmx      = maximum number of l,n,m
    !o zpsi2b     =  the matrix elements
    !r coskt,sinkt= exp(ik.T)
    !h           Ferdi 92.03.17       : original for ASA
    !h           takao at Apr2002(v2) : mixed basis version 2
    !h                    Feb2006(v3) : Add case with ncc/=0
    !r
    !-------------------------------------------------------------------------
    integer:: nl,nlmto,nbloch,nclass,natom,nctot,ncc,nt0,ntp0,mdimx,nlnmx,nnc, &
         ib,ias,ics,ia,ic,nc,nv,nc1,iap,icp,i,itp,it,icr,icore(nl*nl*nnc,nclass),ncore(nclass), &
         mdim(nclass),iclass(natom),imdim(natom),iatomp(natom),verbose,ixx
    real(8)::   ppb(nlnmx,nlnmx,mdimx,nclass) 
    complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc), cphikq(nlmto,*), cphik(nlmto,*),phase(natom)
    if(ncc/=0 .AND. ncc/=nctot) call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
    if(sum(ncore(iclass(1:natom)))/= nctot) call rx( "psicb_v3:sum(ncore) wrong")
    zpsi2b = 0d0
    do concurrent(ia = 1:natom)
       ic    = iclass(ia)
       nc    = nlnmc(ic)
       nv    = nlnmv(ic)
       nc1   = nc + 1
       iap   = iatomp(ia)
       icp   = iclass(iap)
       ib    = imdim(iap)-1
       ias   =  sum(nlnmv(iclass(1:ia-1)))+1
       ics   =  sum(ncore(iclass(1:ia-1)))
       do concurrent(i=1:mdim(icp), itp=1:ntp0, it=1:ncore(ic)) !product basis
          zpsi2b(imdim(iap)+i-1,ics+it,ncc+itp) =phase(ia)*dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icore(it,ic),i,icp)))
       enddo
       if(ncc==0) cycle
       do concurrent(i=1:mdim(icp), itp=1:ncore(ic), it=1:nt0)
          zpsi2b(imdim(iap)+i-1,nctot+it,ics+itp) =dconjg(phase(ia))*sum(cphik(ias:ias+nv-1,it)*ppb(nc1:nc+nv,icore(itp,ic),i,icp))
       enddo
    enddo
  end subroutine psicb_v3
  subroutine drvmelp3(q,ntp0,q_rk,nt0,qik,isp,ginv,ngc,ngcmx,ngpmx,nband,itp0, symope, shtv, qlat, qlatinv,qibz,qbz,nqbz,nqibz, &
       rmel, cmel, nbloch,noccx,nctot, &
       zmelt)
    use m_readqg,only: readqg
    use m_readeigen,only:readgeigw
    ! this is for Wanner (readeigW, drvmelp3)
    implicit none
    real(8):: q(3),q_rk(3),qik(3),ginv(3,3)
    integer:: ngp1, ngp2, ngpmx,nqbz,nqibz, ngcmx ,nctot,nband, &
         ntp0,nt0,nbloch,noccx,  itx, ngc,nnum,inum,ig1,ig2,igc, &
         ngvecpB1(3,ngpmx), &
         ngvecpB2(3,ngpmx), &
         ngveccBr(3,ngcmx), itp0(ntp0), &
         nadd(3),isp  !,ngpn(nqbz)
    complex(8),allocatable::  zmelpl(:,:,:),geigq(:,:)
    real(8):: qlat(3,3),shtv(3),qdiff(3),add(3) &
         ,qibz(3,nqibz),qbz(3,nqbz),qlatinv(3,3),symope(3,3) &
         ,rmel(nbloch,noccx,ntp0) &
         ,cmel(nbloch,noccx,ntp0) &
         ,pi=3.1415926535897932D0
    complex(8):: geig1(ngpmx,nband),geig2(ngpmx,nband) &
         ,zmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
    real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
    integer::verbose
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
    if (nt0 /= ntp0) stop 'drvmelp3: nt0 /= ntp0'
    call melpln2t( ngp1, ngvecpB1 &
         ,   ngp2, ngvecpB2 &
         ,   ngc,  nadd, &
         geig1(1:ngp1,itp0(1:ntp0)), ntp0, & ! &  q1=q    ---> iq ngp1 1:ntp0 q-point
         geig2(1:ngp2,itp0(1:nt0)), nt0,     &! &  q2=q-rk ---> kp ngp2 1:nt0  occupied
         shtv, matmul(symope,qik),qik, symope, qlat, qt,& ! & qt is dummy...
         zmelpl)
    zmelt=0d0
    zmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = &
         dcmplx(rmel (1:nbloch, 1:nctot+nt0, 1:ntp0),  cmel (1:nbloch, 1:nctot+nt0, 1:ntp0))
    zmelt(nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0) =zmelpl(1:ngc,  1:nt0, 1:ntp0)
    deallocate(zmelpl)
  end subroutine drvmelp3
end module m_zmel
