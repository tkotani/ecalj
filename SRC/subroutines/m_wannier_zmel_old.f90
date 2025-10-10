module m_zmel_old !for wannier part This will be removed soon. 
  use m_genallcf_v3,only: natom,nspin,nl,nn,nnv,nnc,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw, nband
  use m_genallcf_v3,only: alat,esmr,nlnmv,nlnmc,icore,ncore,plat,pos,z,ecore,mnl=>nlnm,nl,nn,nlnmx,il,in,im
  use m_hamindex,only: ngrp, symgg=>symops,invg=>invgx
  use m_rdpp,only: Rdpp, nxx,lx,nx,mdimx,nbloch,cgr,ppbrd,nblocha,done_rdpp
  use m_readeigen,only: Readcphif 
  use m_read_bzdata,only: nqbz,nqibz,  qlat,ginv,qbz,qibz,wbz, done_read_bzdata
!  use m_readhbe,only: nband
  use m_itq,only: itq,ntq
  use m_readQG,only: ngpmx,ngcmx,Readqg
  use m_hamindex0,only: Readhamindex0
  use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
  ! OUTPUT: zmel(nbb,nmtot, nqtot) ,nbb:mixproductbasis, nmtot:middlestate, nqtot:endstate
  complex(8),allocatable,protected,public :: zmel(:,:,:)
!  real(8),protected,public:: qm0(3) !for zmel0
  integer,protected:: nbb           !1st dimension of zmel. MPB
  public:: drvmelp3, ppbafp_v2
  
  private
  integer,allocatable,private :: miat(:,:)
  real(8),allocatable,private :: tiat(:,:,:),shtvg(:,:)
  real(8),allocatable,private :: ppbir(:,:,:,:,:,:)
  complex(8),allocatable,private :: ppovlz(:,:)
  real(8),private:: qlatinv(3,3),q_bk(3)=1d10,qk_bk(3)=1d10
  logical,private:: init=.true.
  complex(8),allocatable,private :: cphiq(:,:), cphim(:,:),cphitemp(:,:)
  integer:: nkmin, nkqmin, isp_k, isp_kq,nmtot,nqtot,ispq_bk,ispm_bk
  logical:: debug=.true. !,zzmel0=.false.
  integer:: nbbx=0
  real(8), parameter :: pi=3.1415926535897932D0
contains
  subroutine ppbafp_v2 (irot,ng,isp, mdimx,lx,nx,nxx, &! & Bloch wave
    cgr,lmxax,    & !rotated CG
    ppbrd,        & !radial integrals
    ppb)
    ! <Phi(RLn) Phi(RL'n') B(R,i)>. Atomic part of MPB
    !  n differenciates core phi phidot localOrbital.
    !  in,il,im      = index for n,l,m s. indxlnm.f
    !   ppb            = <Phi(RLn) Phi(RL'n') B(R,i)>
    implicit none
    integer,intent(in) :: irot,ng,isp!,nspin,natom,nlnmx,mdimx
    integer,intent(in) :: lx(natom),nx(0: 2*(nl-1),natom)
    integer,intent(in) :: lmxax,nxx,mdimx
    real(8), intent(in) :: cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng)
    real(8), intent(in) :: ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nxx,natom*nspin)
    real(8), intent(out) :: ppb(nlnmx,nlnmx,mdimx,natom)
    integer :: ic, i,lb,nb,mb,lmb,i1,ibas,i2, np,lp,mp,lmp,n,l,m,lm
    do concurrent (ic=1:natom)
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
  subroutine melpln2t( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd, &
       geigq1, ntp0, &!  q1=q    ---> iq 1:ntp0 q-point
       geigq2, nt0,  &!  q2=q-rk ---> kp 1:nt0  occupied
       shtv,qq, qi, symope, qlat, qt, &
       zmelp)
    use m_read_ppovl,only: getppx2,igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvggg,nvgcgp2,ngvecc
    use m_read_ppovl,only: nggg,ngcgp,ngcread, ggg, ngc2,ngvecc2 !,ppovlinv
    implicit none
    integer:: itp,igc
    integer,intent(in) :: ngp1, ngvecp1(3,ngp1), ngp2
    integer,intent(in) :: ngvecp2(3,ngp2), ngc,nadd(3),ntp0,nt0
    complex(8),intent(in) :: geigq1(ngp1,ntp0),geigq2(ngp2,nt0)
    real(8),intent(in) :: shtv(3),qq(3),qi(3), symope(3,3),qlat(3,3)
    real(8),intent(in) :: qt(3)
    complex(8),intent(out) :: zmelp(ngc,nt0,ntp0)
    complex(8),parameter :: img=(0d0,1d0)
    real(8), parameter :: pi=3.1415926535897932D0
    integer :: nn(1:3)
    integer::igcgp2,iggg,igp1,igp2
!    integer,allocatable:: ngveccR(:,:)
    complex(8),allocatable::ggitp(:,:)
    integer:: ngcgp2,ngcs(1),ngveccR(1:3,1:ngc)
    complex(8)::zdotc
    complex(8),allocatable:: z2(:,:)
    logical:: debug=.false.
    integer:: verbose!,nxi,nxe,nyi,nye,nzi,nze
    complex(8):: zmelp0(ngc,nt0,ntp0),ggitp_(ngp2,ngc),phase(ngc)!,ggitp(ntp0,ngcgp)
    if(ngc==0) return
    if(verbose()>=90) debug= .TRUE. 
    call getppx2(qlat,qi) ! read and allocate PPOVL*
    if(ngc/=ngcread)call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')!  write(6,*)qi,ngcread,ngc
    ngcs(1) = ngc
    call rotgvec(symope, 1, ngc, ngcs, qlat, ngvecc, ngveccR)
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
    ! zmelp <=  \sum_G2 ggitp(Gc+G2) geigqg2(G2)) !! note \bfr'= g (\bfr) +\delta_g  (\bfr= {\bf r})
    ! mapping of function g[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
    ! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
    zmelp0=0d0
    phase(:)=[(exp( img*2d0*pi*sum((qq+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]
    do concurrent (itp= 1:ntp0)
       do concurrent (igc=1:ngc,igp2=1:ngp2)
          nn = ngveccR(:,igc) + ngvecp2(:,igp2)
          igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
          ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
       enddo
       zmelp0(:,:,itp)= matmul(transpose(ggitp_),geigq2)
    enddo
    deallocate(ggitp)
    zmelp=zmelp0 !    call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
    if(verbose()>=100) write(6,*)' melpln2t: end'
  end subroutine melpln2t
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
    zmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = dcmplx(rmel (1:nbloch, 1:nctot+nt0, 1:ntp0),  cmel (1:nbloch, 1:nctot+nt0, 1:ntp0))
    zmelt(nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0) = zmelpl(1:ngc,  1:nt0, 1:ntp0)
    deallocate(zmelpl)
  end subroutine drvmelp3
end module m_zmel_old

