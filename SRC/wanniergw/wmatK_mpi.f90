module m_zmel_old !for wannier part This will be removed soon. 
  use m_genallcf_v3,only: nclass,natom,nspin,nl,nn,nnv,nnc, nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw
  use m_genallcf_v3,only: alat,delta,deltaw,esmr,iclass,nlnmv,nlnmc,icore,ncore,plat,pos,z,ecore,mnl=>nlnm,nl,nn,nlnmx,il,in,im
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
  complex(8),allocatable,protected,public :: zmel(:,:,:)!,zmel0(:,:,:) !for  Get_zmel
  !                                                      zmel0 is just one another zmel which can be contained in m_zmel
  real(8),protected,public:: qm0(3) !for zmel0
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
  logical:: debug=.false. !,zzmel0=.false.
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
  subroutine melpln2t( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd, &
       geigq1, ntp0, &!  q1=q    ---> iq 1:ntp0 q-point
       geigq2, nt0,  &!  q2=q-rk ---> kp 1:nt0  occupied
       shtv,qq, qi, symope, qlat, qt, &
       zmelp)
    use m_read_ppovl,only: getppx2,igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvggg,nvgcgp2,ngvecc
    use m_read_ppovl,only: nggg,ngcgp,ngcread, ggg,ppovlinv, ngc2,ngvecc2
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
    call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
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


subroutine wmatqk_mpi(kount,irot,ef,ef2,esmr,esmr2,tr, &
     ! akao to fit to new hvccfp0.
     ! m, 070501
     !     i          iatomp,rsite,nsp,isp, !ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     iatomp, &
     !     i          rws,irws,nrws,
     rws1,rws2,nrws1,nrws2,nrws, &
     nsp,isp, &! & ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     ifrcw,ifrcwi, &
     qbas,ginv, &
     qibz,qbz,wk,nstbz,wik,nstar,irk, & ! & koun,,iindxk
  
     iclass,mdim,nlnmv,nlnmc, &
     icore,ncore,imdim, &
     ppb, &
     freq_r,freqx,wx,expa,ua,dw,& ! & deltaw,freq
     ecore, &
       
     nlmto,nqibz,nqbz,nctot, &
     !     i          index_qbz, n_index_qbz,  !jan2004
     nl,nnc,nclass,natom, &
     nlnmx,mdimx,nbloch,ngrp,nw_i,nw,nrw,niw,niwx,nq, &
       
     !     &     nblochpmx ,ngpn,ngcni,ngpmx,ngcmx,geigB,ngvecpB,ngveccBr,
     nblochpmx ,ngpmx,ngcmx, &! & ngveccBr,!Jan2004
     wgt0,wqt,nq0i,q0i,symope,alat, shtv,nband, ifvcfpout, &
       !     &     shtw,
     exchange, &! & tote,screen,cohtest, ifexsp,
  ! etra
  ! etra     &     wtet,wtetef,
  ! etra    &     ntqx,ibzx,tetraex,

  !     i omega,iwini,iwend,
  !     i     nbmx,ebmx, !takao 18June2003
     pomatr, qrr,nnr,nor,nnmx,nomx,nkpo,& ! & oct2005 for pomat
     nwf, &
     rw_w,cw_w,rw_iw,cw_iw)
  use m_zmel_old,only: drvmelp3


  ! 2006 May Takashi Miyake, updated for new fpgw
  ! 2004 Sep Takashi Miyake, off-site W
  ! 2004 Jul Takashi Miyake,
  ! 2004 Apr Takashi Miyake, from sxcf_fal2.f

  ! 2001 Sep. esec=omega(itp,iw). Genral iw mode for exchange =F

  ! 2000 takao kotani. This sxcf is modified from sec.f F.Aryasetiawan.

  !  exchange=T : Calculate the exchange self-energy
  !  exchange=F : Calculate correlated part of the self-energy

  !---- correlation case documents by ferdi.Aryasetiawan.  -----------------
  ! 92.02.24
  ! 93.10.18 from sec.f modified to take into account equivalent atoms

  ! the screened coulomb potential
  ! Wc(r,r';w)  = W(r,r';w) - v(|r-r'|)
  !             = < [r1,r2] v(|r-r1|) X(r1,r2;w) v(|r2-r'|) >
  ! W(r,r';w)   = < [r''] ei(r,r'';w) v(|r''-r'| >
  ! ei          = e^(-1), inverse dielectric matrix
  !             = 1 + vX
  ! e           = 1 - vX0 in RPA

  ! expand Wc(r,r';w) in optimal product basis B
  ! Wc(r,r';w)  = S[k=FBZ] S[i,j=1,nbloch]
  !               B(k,i,r) Wc(k,w)(i,j) B(k,j,r')^*
  ! Wc(k,w)(i,j) are  the matrix elements of Wc in B

  ! q       = q-vector in SEc(q,t)
  ! itq     = states t at q
  ! ntq     = no. states t
  ! eq      = eigenvalues at q
  ! ef      = fermi level in Rydberg
  ! tr      = translational vectors in rot*R = R' + T
  ! iatomp(R) = R'
  ! ifrw,ifcw,ifrwi,ifcwi
  !   = direct access unit files for Re and Im coulomb matrix
  !     along real and imaginary axis
  ! ifrb,ifcb,ifrhb,ifchb
  !         = direct access unit files for Re and Im b,hb
  ! qbas    = base reciprocal lattice vectors
  ! ginv    = inverse of qbas s. indxrk.f
  ! xxxx ippb,ipdb,idpb,iddb = pointers to work array w for
  !  ppb     = <phi(RLn) phi(RL'n') B(R,i)>
  !  pdb     = <phi(RLn) phidot(RL'n') B(R,i)>
  !  dpb     = <phidot(RLn) phi(RL'n') B(R,i)>
  !  ddb     = <phidot(RLn) phidot(RL'n') B(R,i)>
  ! freq    = frequencies along real axis
  ! freqx   = gaussian frequencies x between (0,1)
  ! freqw   = (1-freqx)/freqx
  ! wx      = weights at gaussian points x between (0,1)
  ! ua      = constant in exp(-ua^2 w'^2) s. wint.f
  ! expa    = exp(-ua^2 w'^2) s. wint.f
  ! dw      = frequency mesh along real axis
  ! deltaw  = energy mesh in SEc(qt,w) ---Not used now
  ! iclass  = given an atom, tells the class
  ! wk      = weight for each k-point in the FBZ
  ! indexk  = k-point index
  ! qbz     = k-points in the 1st BZ
  ! nstar   = no. stars for each k
  ! irk(k,R) = gives index in the FBZ with k{IBZ, R=rotation
  ! mdim    = dimension of B(R,i) for each atom R
  ! work arrays:
  ! rbq,cbq     = real and imaginary part of b(q)
  ! rhbq,chbq   = real and imaginary part of hb(q)
  ! rbkq,cbkq   = real and imaginary part of b(q-k)
  ! rhbkq,chbkq = real and imaginary part of hb(q-k)
  !   b is the eigenvector of the LMTO-Hamiltonian
  ! ekq     = eigenvalues at q-k
  ! rmel,cmel = real and imaginary part of
  !             <psi(q,t') | psi(q-k,t) B(k,R,i)>
  ! wr1 ... = work arrays
  ! dimensions:
  ! nqibz   = number of k-points in the irreducible BZ
  ! n1,n2,n3= divisions along base reciprocal lattice vectors
  ! natom   = number of atoms
  ! nctot   = no. allowed core states
  ! nbloch  = total number of Bloch basis functions
  ! nlnmx   = maximum number of l,n,m
  ! nlmto   = total number of LMTO basis functions
  ! ngrp    = no. group elements (rotation matrices)
  ! niw     = no. frequencies along the imaginary axis
  ! nw      = no. frequencies along the real axis
  ! niwx    = max(niw,nw)

  !----------------------------------------------------------------------
  use m_readqg,only: readqg0
  use m_readeigen,only:readcphiw
  use m_keyvalue,only: getkeyvalue
  use m_read_bzdata,only: wklm

  ! RS: modules for MPI
  use rsmpi
  use rsmpi_qkgroup
  use rsmpi_rotkindex
  implicit none
  integer(4) :: ntq, natom,nqbz,nqibz,ngrp,nq,nw_i,nw,niw, &
       nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,nstatex, &
       irot,  iqisp,ikpisp,isp,nsp,  nlnmx, iq, ip, it,itp, it2, itp2,  iiclass,mdim(*), &
       ifrcw,ifrcwi, ifvcfpout,ndummy1,ndummy2,kx,kr,kr2,kr3,ngc,ngb,nbloch, &
       kp,nt0,nocc, nt0p,nt0m,irkp,i,nt0org,nmax,nt,ntp0, &
       nbmax,nclass,nl,nnc, nblochpmx,ix,nx,iw,iwp,ixs,ixsmx, &
       mdimx, nwx,niwx, iatomp(natom), &
       nstar(nqibz),irk(nqibz,ngrp),kount(nqibz),nwf !,iclose
  real(8) :: q(3),qbas(3*3),ginv(3*3),tr(3,natom), &
       wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz), &
       freqx(niw),wx(niw),expa(niw), &
       eq(nband), &
       !     &   ekq(nband), ekc(nctot+nband),
       tpi,ef,ef2,esmr,esmr2,efp,efm,wtx,wfac,wfacx,we,esmrx,ua, &
       dw,wtt,wexx,www,exx,exxq,weight
  !      complex(8) :: zsec(-1:1,ntq,nq)
  !      real(8)    ::  shtw
  !                       ! This shft is  to avoid some artificial resonance effects.
  !                       ! shtw can be zero for esmr/=0 given by takao.

  integer(4):: ngpmx, ngcmx,  igc, nadd(3)
  real(8) :: wgt0(nq0i,ngrp),wqt(nq0i),qk(3), qbasinv(3,3),qdiff(3),add(3),symope(3,3), &
       qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot), &
       ppb(1) !pdb(1),dpb(1),ddb(1)
  real(8),allocatable:: rmelt(:,:,:),cmelt(:,:,:), rmelt2(:,:,:),cmelt2(:,:,:), &
       rmelt3(:,:,:,:),cmelt3(:,:,:,:)
  complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:), &
       zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:), &
       zmelt(:,:,:),zmelc(:,:,:,:),zmelcc(:,:,:,:)
  ! for exchange --------------------
  logical :: exchange,screen,cohtest,tote
  real(8),allocatable:: &
       w1p(:,:,:),w2p(:,:,:)
  complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)

  !- debug write ---------------------
  logical :: debug=.false.

  integer(4) :: ibl,iii,ivsumxxx,ifexsp ,iopen
  integer(4),save::ifzwz=-999

  integer(4) :: iwini, iwend, ia
  !      real(8)    :: esec, omega(ntq, iwini:iwend)
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       cw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       rw_iw(nwf,nwf,nwf,nwf,nrws,niw), &
       cw_iw(nwf,nwf,nwf,nwf,nrws,niw)
  complex(8),allocatable:: expikt(:)
  complex(8):: img=(0d0,1d0)
  ! akao
  complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband) &
       , cphiqtmp(nlmto,nband)

  ! cccccccccccccccccccccccccccccccccccccccccccccc faleev 2002
  integer(4) :: nt_max, igb1,igb2,iigb, nw_w
  complex(8),allocatable:: zmel1(:)
  complex(8), allocatable :: zw_(:,:) !,zzmel(:,:)
  complex(8), allocatable :: zwz2(:,:),zw2(:,:,:,:) !0 variant
  complex(8) ::  zz2 ,zwz3(3)
  real(8) :: dd,omg_c,dw2
  real(8) :: freq_r(nw_i:nw)
  complex(8), allocatable :: zw3(:,:,:)


  real(8)::weavx,wfaccut=1d-10

  logical :: GaussSmear
  real(8) :: ddw !ebmx,
  integer(4):: nbmxe,nstatetot !nbmx,

  !      integer(4):: n_index_qbz
  !      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

  integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

  integer(4)::verbose,nstbz(nqbz),iqini,iqend !bzcase,
  real(8):: wgtq0p

  integer(4):: iqindx,nrec,kxx
  real(8)::quu(3),qibz_k(3),qbz_kr(3)
  logical :: onlyQ0P, onlyimagaxis ,noq0p !,noq0p,test_omitq0p,

  logical ::zwz3mode
  !      logical ::testimx=.false.

  real(8):: ua_,expa_(niw),ua2,freqw,freqw1,ratio,ua2_(niw)
  logical :: ua_auto
  integer(4):: icc=0
  real(8),allocatable:: uaa(:,:)

  !      logical ::testimx=.false.
  ! ccc zvz test cccccccccccccccccccccccccc
  integer(4):: ngbx
  !      complex(8):: vcoul(ngbx,ngbx)
  complex(8),allocatable:: vzz(:,:,:),aaa(:)
  complex(8):: zvz,zvz1
  integer(4):: ib1,ib2,ifix
  ! ccccccccccccccccccccccccccccccccc
  integer(4) ::nbcut,nbcutc
  logical ::iww2=.true., oncew


  !...
  logical::smbasis
  integer(4):: nn,no,ifpomat,iclose,isx,iqx
  complex(8),allocatable:: pomat(:,:)
  real(8):: q_r(3)
  integer(4):: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo)
  complex(8):: pomatr(nnmx,nomx,nkpo)
  real(8):: qrr(3,nkpo)

  real(8):: elxx,ehxx,ekxx,efxx
  integer(4):: ixsmin,iwm,iir,nwxi
  real(8)   :: fffr(3)
  complex(8):: zwzz(3)

  ! m
  integer(4) :: nqbz2,nwf2,iko_ix,iko_fx,iqtmp,ifmlw,nko,iqk &
       ,ifi,in1,in2,imp,ilp,ii,jj,nrws,nrws1,nrws2 &
       ,ir1,ir2,ir3,ir,nrw
  real(8) :: norm2,qtmp(3),rws1(3,nrws1),rws2(3,nrws2),tmp
  complex(8) :: ztmp,expiqR1(nrws1),expiqR2
  complex(8),allocatable :: cnk(:,:,:),zmel2(:,:,:),zmel3(:,:,:)
  integer(4) :: itq(nwf)
  complex(8) :: weightc(nrws1),zmeltt1

  complex(8),allocatable:: ppovl(:,:),ppovlz(:,:),zcousq(:,:),zmeltt(:,:,:)
  real(8),allocatable::vcoud(:),vcousq(:)
  real(8),parameter:: pi=4d0*atan(1d0),fpi=4d0*4d0*atan(1d0)

  integer:: ifvcoud,ivc,ngb0
  real(8)::qvv(3),vc
  logical :: newansisoW
  character(10):: i2char

  !      real(8),allocatable:: wklm(:),dmlx(:,:),epinvq0i(:,:)
  !      integer:: ifidmlx,lxklm,nq0ix
  complex(8)::w3p

  integer:: mrecl,nprecx,nwordr,il

  !--------------------------------------------------------------------
  ! RS: MPI parameters
  integer :: kx_local
  !--------------------------------------------------------------------
  debug=.false.
  if(verbose()>=90) debug= .TRUE. 

  !     ! === readin Vcoud and EPSwklm for newansisoW()=T ===
  !$$$      if(newansisoW()) then
  !$$$         ifidmlx = iopen('EPSwklm',0,0,0)
  !$$$         read(ifidmlx) nq0ix,lxklm
  !$$$         if(nq0i/=nq0ix) then
  !$$$            write(6,*)'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
  !$$$            call rx( 'nq0i from EPSwklm /= nq0i')
  !$$$         endif
  !$$$         allocate( dmlx(nq0i,9))
  !$$$         allocate( epinvq0i(nq0i,nq0i) )
  !$$$         allocate( wklm((lxklm+1)**2))
  !$$$         read(ifidmlx) dmlx, epinvq0i
  !$$$         read(ifidmlx) wklm
  !$$$c         do il=1,(lxklm+1)**2
  !$$$c         write(6,"('EPSwklm=', i5,3f12.5)") il,wklm(il)
  !$$$c         enddo
  !$$$         ifidmlx = iclose('EPSwklm')
  !$$$      endif


  ! oct2005
  call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
  nbcutc=nctot+nbcut

  tpi         = 8d0*datan(1.d0)
  !      iq         = idxk (q,qbz,nqbz) ! index for q
  !      write(6,"(' iq q  =',i4,3f8.4)")iq,q
  ! cc      iq          = idxk (q,qbze,nqbze) ! index for q
  !      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
  nlmtobnd    = nlmto*nband
  nstatetot      = nctot + nband
  call minv33(qbas,qbasinv)

  ! work arrays for psi2br.f
  if(debug) write(6,*) ' sxcf: 1'
  allocate(expikt(natom))

  !      if(bzcase()==1) then
  if(abs(sum(qibz(:,1)**2))/=0d0) call RSMPI_Stop( ' sxcf assumes 1st qibz/=0 ')
  if(abs(sum( qbz(:,1)**2))/=0d0) call RSMPI_Stop( ' sxcf assumes 1st qbz /=0 ')
  !      endif

  do it = 1,nwf
     itq(it) = it
  enddo

  ! m debug
  !        write(*,*)'isp,',isp,rw_w(1,1,1,1,1,0),rw_w(2,2,2,2,1,0)

  !-----
  if(exchange .AND. ( .NOT. newansisoW())) then
     rewind  ifvcfpout
     read(ifvcfpout) ndummy1, ndummy2
  endif

  !===============================
  ! loop over irreducible k-points
  !===============================
  ! ccccccccccccccccccccccccccccccc
  !      iii = ivsumxxx(irk,nqibz*ngrp)
  !      write(6,*)' sxcf:sum non-zero irk=',iii
  !      stop "sss"

  ! ccccccccccccccccccccccccccccccc

  !      if(bzcase()==1) then
  kx = 1  ! qibz(:,1)=0 contribution for kcount
  if(irk(kx,irot)/=0) kount(kx)= kount(kx) + 1
  !      kount(kx)= kount(kx) + 1
  !      endif

  ! --- main loop start
  iqini=2
  !      if(bzcase()==2) iqini=1
  iqend=nqibz+nq0i

  if(newansisoW()) then       !takao2012apr
     iqini=1
     iqend=nqibz            !no sum for offset-Gamma points.
     !     if(exchange) iqend=nqibz
  endif
  ! ccccccccccccccccccccccccc
  !      iqini=2
  ! ccccccccccccccccc

  ! RS: ifile_rsmpi is defined in gwsrc/RSMPI_mod.F
  write(ifile_rsmpi,*) "RS: irot = ", irot
  write(ifile_rsmpi,*) "RS: main loop start"

  ! cccccccccccccccccccccccccccc
  call getkeyvalue("GWinput","TestOnlyQ0P",onlyq0p,default=.false.)
  call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
  if ( .NOT. noq0p) &
       call getkeyvalue("GWinput","NoQ0P",noq0p,default= .FALSE. )
  if (noq0p)write(*,*)'noq0p mode'
  if(noq0p) iqend=nqibz
  !      iqend=nqibz
  !      if(test_omitq0p()) then
  !        iqend=nqibz
  !        write(6,*)'iqend=',iqend
  !      endif
  ! cccccccccccccccccccccccccc
  do 1100 kx_local = 1,nk_local_rotk(irot)
     kx = ik_index_rotk(irot,kx_local)
     !      do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
     ! debug:
     !      do 1100 kx = iqini,iqini !kx=1 corresponds to q=0 is omitted.
     !        write (6,"(i3,'  out of ',i3,$)") kx,iqend
     write(ifile_rsmpi,*) ' wmat_MPI: goto loop kx=',kx
     if(debug)  write(6,*) ' sxcf: goto loop kx=',kx

     !        write(*,'("1  begin k-cycle",$)')
     !         call cputid(0)
     !          write(*,*)'kx, ip, irot=',kx, ip,irot

     if( kx <= nqibz ) then
        !          k  = kx
        kr = irk(kx,irot) ! index for rotated k in the FBZ
        qibz_k= qibz(:,kx)
        !          qbz_kr= qbz (:,kr)
        if(kr/=0) qbz_kr= qbz (:,kr) !feb2006
     else
        !          k = 1  ! corresponds to q=0
        !          kr= 1  ! corresponds to q=0
        !          k = iqindx((/0d0,0d0,0d0/), ginv, qibz,nqibz)
        !          kr= iqindx((/0d0,0d0,0d0/), ginv, qbz,  nqbz)
        kr=-99999 !for sanity check
        qibz_k= 0d0
        qbz_kr= 0d0
     endif
     !        ngc = ngcni(k)  ! k-points in IBZ
     !        write(6,*) ' k ngc=',k,ngc
     !        ngb = nbloch + ngcni(k)

     call readqg0('QGcou',qibz_k,  quu,ngc)
     !        ngc = ngcni(k)  ! k-points in IBZ
     ngb = nbloch + ngc

     !! ===Readin diagonalized Coulomb interaction===
     !! note sep102012takao
     !!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
     !!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
     !!  Vould file is written in hvccfp0.m.F.
     !! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
     !! These can be also separeted into WVR.ID and WVI.ID files.
     if(newansisoW()) then
        qxx=qibz_k
        !           if(kx<=nqibz) qxx=qibz_k
        !           if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
        ifvcoud = iopen('Vcoud.'//i2char(kx),0,0,0)
        do
           read(ifvcoud) ngb0
           read(ifvcoud) qvv
           !              write(6,"('readin qvv ngb0=',3f9.4,i5)")qvv,ngb0
           !              write(6,"('readin qxx ngb0=',3f9.4,i5)")qxx
           if(allocated(vcoud)) deallocate(vcoud)
           allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
           read(ifvcoud) vcoud
           !$$$  cccccccccccccccccccccccccccccccccccccccccccccc
           !$$$  if(sum(abs(qxx(1:2)))<1d-4) then
           !$$$  do irot2 = 1,ngrp
           !$$$  kr = irkip(kx,irot2,ip) ! index for rotated kr in the FBZ
           !$$$  if(kr==0) cycle ! next irot
           !$$$  vcoud(1) =vcoud(1)*wqfac(kr)
           !$$$  exit
           !$$$  enddo
           !$$$  write(6,*)'wwwwwwwwwww why here ?wwwwwwwwww'
           !$$$  endif
           !$$$  cccccccccccccccccccccccccccccccccccccccccccccc
           read(ifvcoud) zcousq
           if(sum(abs(qvv-qxx))<1d-6) goto 1133
        enddo
        if(sum(abs(qvv-qxx))>1d-6) then
           write(6,*)'qvv =',qvv
           write(6,*)'qxx=',qxx,kx
           call rx( 'wmatK: qvv/=qibz(:,kx) hvcc is not consistent')
        endif
1133    continue
        ifvcoud = iclose('Vcoud.'//i2char(kx))
        if( ngb0/=ngb ) then !sanity check
           write(6,*)' qxx ngb0 ngb=',qxx,ngb0,ngb
           call rx( 'hsfp0.m.f:ngb0/=ngb')
        endif
        !$$$  if(sum(abs(qibz_k))<1d-6) then
        !$$$  idummy  = iclose('Vcoud') !close and open again. This is because first data in Voud is for q=0
        !$$$  ifvcoud = iopen('Vcoud',0,0,0)
        !$$$  endif
        !$$$  read(ifvcoud) ngb0
        !$$$  if( ngb0/=ngb ) then
        !$$$  write(6,*)' qibz_k=',qibz_k,ngb0,ngb
        !$$$  write(6,*)' qibz_k=',qibz_k
        !$$$  stop 'hsfp0.m.f:ngb0/=ngb'
        !$$$  endif
        !$$$  read(ifvcoud) qvv
        !$$$  if(sum(abs(qvv-qibz_k))>1d-6) then
        !$$$  write(6,*)'qvv =',qvv
        !$$$  write(6,*)'qibz_k=',qibz_k,kx
        !$$$  stop 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent'
        !$$$  endif
        !$$$  if(allocated(zcousq)) deallocate(zcousq,vcousq,vcoud)
        !$$$  allocate( zcousq(ngb0,ngb0),vcousq(ngb0),vcoud(ngb0) )
        !$$$  read(ifvcoud) vcoud
        !$$$  read(ifvcoud) zcousq
        !$$$  vcousq=sqrt(vcoud) !

        !! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
        !! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi_q,itp |phi_q-rk,it B_rk,I> O^-1_IJ
        !! ppovlz= O Zcousq
        !! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
        if(allocated(ppovlz)) deallocate(ppovlz)
        allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
        call readppovl0(qibz_k,ngc,ppovl)
        ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
        ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
        deallocate(zcousq,ppovl)
     endif

     if( .NOT. newansisoW()) then
        if(exchange) then
           ! RS: set pointer to the right place
           call set_vcoul_rsmpi(ifvcfpout,kx-iqini+1)
           !     allocate(vcoul(ngb,ngb))
           !     read(ifvcfpout) vcoul(1:ngb,1:ngb)
           read(ifvcfpout) nn !oct2005
           allocate(vcoul(nn,nn))
           read(ifvcfpout) vcoul(1:nn,1:nn)
        endif
        !     ! weight check for cycle or not.
        if(kx <= nqibz ) then
           if (kr == 0)    then
              !     stop 'wmat: kr=0'
              if(exchange) deallocate(vcoul)
              cycle
           endif
           kount(kx)= kount(kx) + 1 ! count the no. times k
           ! appears in the 1st BZ
           ! cccccccccccccccccccccccccccccccccccccccccccccccc
           !     write(6,*)' irot,ip, k, kount in  =',irot, ip, k, kount(k,ip)
           !     deallocate(vcoul)
           !     cycle
           !     write(6,*)' kount out =',kount(k)
           ! ccccccccccccccccccccccccccccccccccccccccccccccccc
           !     if (kount(kx) > nstar(kx)) stop 'sexc: too many stars'
           if (kount(kx) > nstar(kx)) stop 'wmat: kount > nkstar'
           !     if (kount(kx) > 1) stop 'wmat: kount > 1'
        else
           if( wgt0(kx-nqibz,irot) == 0d0 ) then
              if(exchange) deallocate(vcoul)
              cycle
           endif
        endif
     else
        if (kr == 0) cycle
     endif
     !---test
     if(OnlyQ0P .AND. kx<=nqibz) then
        if(exchange) deallocate(vcoul)
        cycle
     endif
     !! phase factor for off-site W
     do ir1=1,nrws1
        expiqR1(ir1) = exp(-img*tpi* sum(qbz_kr(:)*rws1(:,ir1)))
     enddo


     !! ===================================================================
     allocate( rmelt3(ngb,nwf,nwf,nrws2),cmelt3(ngb,nwf,nwf,nrws2))
     rmelt3 = 0d0
     cmelt3 = 0d0
     !! loop over FBZ
     do iq = 1,nqbz
        q(:) = qbz(:,iq)
        call readcphiW (qbz(:,iq), nlmto,isp, quu, cphiq)
        qk =  q - qbz_kr          ! qbz(:,kr)
        call  readcphiW(qk, nlmto,isp, quu, cphikq)
        do ia = 1,natom
           expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) !  write(6,'(" phase ",i3,2d12.4)')ia,expikt(ia)
        end do
        if(debug) write(6,*) ' sxcf: tr=',tr
        if(debug) write(6,*) ' sxcf: goto psicb2'
        nbmax = nwf
        nt   = nctot + nbmax      ! = nstate for the case of correlation
        ntp0 = nwf
        allocate( zzmel(nbloch,nt,ntp0)) ! rk,ibloch  q-rk,it  q,itp
        zzmel = 0d0
        call psi2b_v2 (nbmax, ntp0, iclass, &
             dreal(expikt(1:natom)),dimag(expikt(1:natom)), &
             cphikq,      &        ! & rbkq,cbkq,rhbkq,chbkq, !  q-rk nko
             cphiq,            &   ! & rbq,cbq,rhbq,chbq,     !  q    nko
             ppb,              &   ! & pdb,dpb,ddb,
             nlnmv,nlnmc,mdim,nctot, &
             imdim,iatomp, &
             mdimx,nlmto,nbloch,nlnmx, nband, nt,ntp0, &
             natom,nclass, &
             zzmel)               ! rmel,cmel)
        if(debug) write(6,"('sum of zmel abszmel=',4d23.16)") &
             sum(zzmel),sum(abs(zzmel) )
        allocate( zmelt(ngb, nctot+nbmax, ntp0) ) ! rk,ibloch  q-rk,it  q,itp
        if(debug) write(6,*) ' sxcf_fal2: goto drvmelp xxxxx1',ngb,nctot,nbmax,ntp0
        call drvmelp3( q,   ntp0, &! &  q in FBZ  q,itp
             q-qbz_kr, nbmax,  &! &  q-rk,it
             qibz_k,           &! &  k in IBZ for e-product basis
             isp,ginv, &
             ngc,ngcmx,ngpmx,nband,itq, &
             symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz, &
             dreal(zzmel), dimag(zzmel), nbloch, nt,nctot, &
             zmelt)
        deallocate(zzmel) !rmel,cmel)
        if(debug) write(6,*) ' sxcf: goto wtt'
        if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)") sum(zmelt)
        do ir2=1,nrws2
           expiqR2 = exp( img*tpi* sum(q(:)*rws2(:,ir2)))
           rmelt3(:,:,:,ir2) = rmelt3(:,:,:,ir2) + wk(iq) * dreal(zmelt(:,:,:)*expiqR2)
           cmelt3(:,:,:,ir2) = cmelt3(:,:,:,ir2) + wk(iq) * dimag(zmelt(:,:,:)*expiqR2)
        enddo  ! ir2
        deallocate(zmelt)
     enddo
     !! ===================================================================



     if(kx<= nqibz) then
        wtt = wk(kr)           !         wtx = 1d0
     else
        wtt = wk(1)*wgt0(kx-nqibz,irot) ! wtx = wgt0(kx-nqibz,irot)
        if(abs(wk(1)-1d0/dble(nqbz))>1d-10) stop 'sxcf:wk(1) inconsistent'
     endif
     weight = wtt
     if(debug) then
        write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
     endif
     do ir1=1,nrws1
        weightc(ir1) = weight*expiqR1(ir1)
     enddo

     !--------------------------------------------------------
     ! --- bare Coulomb section ---
     !--------------------------------------------------------

     ! S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
     !                        v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>

     !> z1p(j,t,t') = S[i=1,nbloch] <psi(q,t') | psi(q-rk,t) B(rk,i)> v(k)(i,j)


     !      write(6,*)' vcoulsum=',sum(vcoul)
     !      if(debug) write(6,*)'  sumz=',dcmplx(rmelt,cmelt),sum(vcoul)

     if(exchange) then
        if (debug) write(*,*) 'bare coulomb section begins'
        allocate(zmel1(ngb))
        allocate(zmel(ngb, nwf, nwf))
        if( .NOT. newansisoW()) allocate(vcoult(1:ngb,1:ngb),z1p(ngb,nwf,nwf))
        do ir2=1,nrws2
           !!  (rmelt3,cmelt3)  !rk,ibloch  q-rk,it  q,itp
           zmel  = dcmplx (rmelt3(:,:,:,ir2),cmelt3(:,:,:,ir2)) !<psi_itp|psi_it B>
           if( .NOT. newansisoW() ) then
              vcoult= transpose(vcoul)
              call matm( vcoult, zmel, z1p,  ngb,ngb,nwf*nwf ) !z1p= <psi_itp|psi_it B> * voul
              !            deallocate(vcoult)
              do ir3=1,nrws2
                 do itp2 = 1,nwf
                    do it2  = 1,nwf
                       do it   = 1,nwf
                          do itp  = 1,nwf
                             zmel1(:)=dcmplx(rmelt3(:,it,itp,ir3),-cmelt3(:,it,itp,ir3)) ! <B psi_it|psi_itp>
                             ztmp = sum ( z1p(:,it2,itp2)*zmel1 ) ! <psi_itp2|psi_it2 B>*vcoul*<B psi_it|psi_itp>
                             do ir1=1,nrws1
                                ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
                                rw_w(itp2,it2,it,itp,ir,0) = rw_w(itp2,it2,it,itp,ir,0) &
                                     + real(ztmp*weightc(ir1))
                                cw_w(itp2,it2,it,itp,ir,0) = cw_w(itp2,it2,it,itp,ir,0) &
                                     + imag(ztmp*weightc(ir1))
                             enddo ! ir1
                          enddo
                       enddo
                    enddo
                 enddo
              enddo ! ir3
           else
              ! based on E_I basis. See Christoph's paper
              allocate(zmeltt(nwf,nwf,ngb))
              do itp= 1,nwf
                 do it = 1,nwf
                    do ivc=1,ngb
                       zmeltt(it,itp,ivc)=sum(zmel(:,it,itp)*ppovlz(:,ivc)) ! <psi_itp|psi_it E_I> (I=ivc)
                    enddo
                 enddo
              enddo
              do ir3=1,nrws2
                 do itp2 = 1,nwf
                    do it2  = 1,nwf
                       do it   = 1,nwf
                          do itp  = 1,nwf
                             zmel1(:)=dcmplx(rmelt3(:,it,itp,ir3),cmelt3(:,it,itp,ir3)) ! <psi_itp|psi_it B_I>
                             w3p=0d0
                             do ivc=1,ngb
                                zmeltt1 =  sum( zmel1(:)*ppovlz(:,ivc) ) !<psi_itp|psi_it E_I>
                                if(ivc==1 .AND. kx==iqini) then
                                   vc= wklm(1)* fpi*sqrt(fpi) /wk(kx) !kx right?
                                else
                                   vc= vcoud(ivc)
                                endif
                                w3p=w3p + zmeltt(it2,itp2,ivc)  *vc* dconjg(zmeltt1)
                                ! <psi_itp2|psi_it2 E_I> *vc* <E_I psi_it|psi_itp>
                             enddo
                             ztmp= w3p
                             do ir1=1,nrws1
                                ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
                                rw_w(itp2,it2,it,itp,ir,0) = rw_w(itp2,it2,it,itp,ir,0) &
                                     + real(ztmp*weightc(ir1))
                                cw_w(itp2,it2,it,itp,ir,0) = cw_w(itp2,it2,it,itp,ir,0) &
                                     + imag(ztmp*weightc(ir1))
                             enddo ! ir1
                          enddo
                       enddo
                    enddo
                 enddo
              enddo ! ir3
              deallocate(zmeltt)
           endif
        enddo ! ir2
        if(allocated(vcoul)) deallocate(vcoul,vcoult)
        if(allocated(z1p))   deallocate(z1p)
        if(allocated(rmelt3)) deallocate(rmelt3,cmelt3,zmel1, zmel)
        if (debug) write(*,*) 'bare coulomb section finished'
        !! --- End of bare-Coulomb section --------------
     else
        !--------------------------------------------------------------------------
        !--- screening effect section----------------------------------------------
        !--------------------------------------------------------------------------
        ! S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
        !                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n') |psi(q,t')>
        ! t-->itp   n -->it
        ! t'-->itp2 n'-->it2
        !--------------------------------------------------------------
        !!--- The matrix elements zmelc.
        !! zmelc  = < E(rk,j) psi(q-rk,it) | psi(q,itp) >
        !! E basis is ghe Christoph's basis diagonalize the Coulomb inteaction
        allocate( zmelc (ngb, nwf, nwf,nrws2),zmelcc (ngb, nwf, nwf,nrws2), &
             zw (nblochpmx,nblochpmx), &
             zw2(nwf,nwf,nwf,nwf) )
        zmelcc = dcmplx (rmelt3,-cmelt3) !zmelcc = < B(rk,j) psi(q-rk,it) | psi(q,itp) >
        do it =1,nwf
           do itp=1,nwf
              do ir2=1,nrws2
                 zmelc(:,it,itp,ir2) =  matmul(zmelcc(:,it,itp,ir2),dconjg(ppovlz(:,:)))
              enddo
           enddo
        enddo
        deallocate(rmelt3,cmelt3,zmelcc)
        if(debug) write(6,*)' end of zmel'
        !====================================================================
        ! Wc(qt,w) along the imaginary axis
        !====================================================================
        !------------------------------------------------
        ! loop over w' = (1-x)/x, frequencies in Wc(k,w')
        ! {x} are gaussian points between (0,1)
        !------------------------------------------------
        nx  = niw
        nprecx=8
        mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
        ifrcwi = iopen('WVI.'//i2char(kx),0,-1,mrecl)
        do ix = 1,nx     ! imaginary frequency w'-loop
           ! ccccccccccccccccccccccccccccccc
           !          nrec=(kx-2)*niw+ix
           !          if(bzcase()==2) nrec= (kx-1)*niw+ix
           nrec=ix
           read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis
           !          read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! Readin W-v on imag axis
           ! ccccccccccccccccccccccccccccccc

           ! zwz= S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
           !                Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
           !        do itp = 1,ntp0
           !        do  it = 1,nstate
           !          zwz(ix,it,itp) = sum(
           !     &   dconjg(zmel(:,it,itp)),matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)) )
           !        enddo
           !        enddo
           do ir3=1,nrws2
              do ir2=1,nrws2
                 call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), zmelc(:,:,:,ir3), &
                      nwf,nwf,ngb, &
                      zw2)
                 do ir1=1,nrws1
                    ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
                    rw_iw(:,:,:,:,ir,ix) = rw_iw(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
                    cw_iw(:,:,:,:,ir,ix) = cw_iw(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
                 enddo ! ir1
              enddo ! ir2
           enddo ! ir3
        enddo
        ifrcwi = iclose('WVI.'//i2char(kx))

        !====================================================================
        ! Wc(qt,w) along the real axis
        !====================================================================
        if(debug) write(6,*)' go to poles'
        nprecx=8
        mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
        ifrcw = iopen('WVR.'//i2char(kx),0,-1,mrecl)
        do      ix = 0,nrw                    ! real frequency w'-loop
           ! ccccccccccccccccccccccccc
           !          nrec=(kx-2)*(nw+1-nw_i)+ ix-nw_i+1
           !          if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
           nrec=ix-nw_i+1
           read(ifrcw,rec=nrec) zw
           ! cccccccccccccccccccccccc
           ! zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
           do ir3=1,nrws2
              do ir2=1,nrws2
                 call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), zmelc(:,:,:,ir3), &
                      nwf,nwf,ngb, &
                      zw2)
                 do ir1=1,nrws1
                    ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
                    rw_w(:,:,:,:,ir,ix)  = rw_w(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
                    cw_w(:,:,:,:,ir,ix)  = cw_w(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
                 enddo ! ir1
              enddo ! ir2
           enddo ! ir3
        enddo
        ifrcw = iclose('WVR.'//i2char(kx))
        deallocate(zmelc,zw,zw2)
        if(debug) write(6,*)' end of screen-if'
        ! end of if (exchange)
     endif
     continue  ! end of k-loop
1100 enddo
  return
  !     end subroutine wmatqk_mpi

  !!-------------------------------------------------------------------------------
contains

  !! psi2b_v2 and psicb_v2 are older versions of psi2b_v3 and psicb_v3 in m_zmel
  subroutine psi2b_v2(nt0,ntp0,iclass,coskt,sinkt, &
       cphik, cphikq,    ppb,  nlnmv,nlnmc,mdim,nctot,imdim,iatomp, &
       mdimx,nlmto,nbloch,nlnmx,noccxv,nt,ntp, &
       natom,nclass, &
       zpsi2b)
    ! originaly 92.03.17 by Ferdi.
    ! takao modified at Apr 2002
    ! calculates <psi(k',t') | psi(k,t) B(R,i)>
    ! for all R
    ! psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
    ! B(R,i)   = Bloch orthonormal product basis for atom R
    ! psi(k,t) is stored after nctot

    ! nt0        = no. t
    ! ntp0       = no. t'
    ! coskt,sinkt= exp(ik.T)
    ! cphik b(k)
    ! cphikq b(k')

    ! ppb        = <phi(RLn) phi(RL'n') B(R,i)>

    ! ddb        = <phidot(RLn) phidot(RL'n') B(R,i)>, s. ppbl.f
    ! nlnmv      = number of l,n,m for valence
    ! nlnmc      = number of n,l,m for core states
    ! mdim       = number of optimal product basis functions for each class
    ! nctot      = total no. allowed core states
    ! nbloch     = total no. optimal product basis
    ! nlnmx      = maximum number of l,n,m
    ! noccxv     = max. no. occupied valence states
    ! nt         = maximum number of occupied states
    ! ntp        = maximum number of unoccupied states

    ! zpsi2b     =  the matrix elements
    implicit real*8(a-h,o-z)
    implicit integer (i-n)
    integer::nt0,ntp0,natom,nlmto,nclass,nlnmc,mdim,nctot
    complex(8):: cphik(nlmto,noccxv),cphikq(nlmto,ntp0) &
         ,zpsi2b(nbloch,nt,ntp),phase
    dimension &
         ppb(nlnmx,nlnmx,mdimx,nclass), &
         nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom), &
         coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)

    integer(4),allocatable::iasx(:)
    integer(4) :: nzppb1,nzppb2
    complex(8),allocatable :: zz(:,:), zppb(:,:)
    complex(8) :: alpha,beta

    ! zppb is used as work array for ppb(:,:,i,ic) and for zpsi2b(ib,:,:).
    nzppb1=max(nt0,nlnmx)
    nzppb2=max(ntp0,nlnmx)
    allocate( zz(nlnmx,ntp) )
    allocate( zppb(nzppb1,nzppb2) )
    beta=0d0  ; alpha=1d0

    ! check dimensions
    if (ntp0 > ntp) call rx( 'psi2bc: ntp exceeded')
    if (mdimx /= maxval(mdim)) call rx( 'psi2bc: wrong mdimx')
    if (nctot+nt0 > nt) call rx( 'psi2bc: nt exceeded')
    if (nt0 > noccxv) call rx( 'psi2bc: noccxv exceeded')
    if ( sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v2: wrong nbloch')
    allocate(iasx(natom))
    ias = 1
    do ia = 1,natom
       iasx(ia) = ias
       ias = ias + nlnmv(iclass(ia))
    enddo
    if(ias-1/=nlmto) call rx( ' psi2b_v2:sum(nlnmv)/= nlmto')
    ! loop over atoms
    do  ia = 1,natom
       ic   = iclass(ia)
       nc   = nlnmc(ic)
       nv   = nlnmv(ic)
       nc1  = nc + 1
       if (nc+ nlnmv(ic) > nlnmx) call rx( 'psi2b_v2: nlnmx exceeded')
       phase= dcmplx(coskt(ia),sinkt(ia))
       ias  = iasx(ia)
       iap  = iatomp(ia)
       icp  = iclass(iap)
       do   i = 1,mdim(icp) ! loop over optimal product basis
          !---------------------------------------------------
          !c sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)>
          !c for a given i, for all L'n' and t'
          !c bkq is complex but < > is real
          !1      do     itp = 1,ntp0
          !1      do      jp = 1,nlnmv(ic)
          !1      zz (jp,itp)=dconjg(
          !1     &    sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,nc+jp,i,icp)) )
          !1      end do
          !1      end do

          !2      zz(1:nv,1:ntp0) =dconjg(
          !2     & matmul(  transpose(ppb(nc1:nc+nv,nc1:nc+nv,i,icp))
          !2     &         ,cphikq(ias:ias+nv-1,1:ntp0)) )

          !3        call dgemm('T','N',nv,ntp0,nv,
          !3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
          !3     &          dreal(cphikq(ias:ias+nv-1,1:ntp0)), nv,
          !3     &          0d0,
          !3     &   rr, nlnmx )
          !3        call dgemm('T','N',nv,ntp0,nv,
          !3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
          !3     &          dimag(cphikq(ias:ias+nv-1,1:ntp0)), nv,
          !3     &          0d0,
          !3     &   cc, nlnmx )
          zppb(1:nv,1:nv) = ppb(nc+1:nc+nv,nc+1:nc+nv,i,icp)
          call zgemm('T','N',nv,ntp0,nv, &
               alpha, zppb,nzppb1, cphikq(ias,1), nlmto,  beta, &
               zz,  nlnmx )
          do itp = 1,ntp0
             do jp = 1,nv
                zz(jp,itp)= dconjg(zz(jp,itp) )
             enddo
          enddo
          !----------------------------------------------------
          ! <psi(k+q,t') | psi(k,t) B(i)>
          !1      do      it = 1,nt0
          !1      do     itp = 1,ntp0
          !1       zpsi2b(ib,nctot+it,itp)=
          !1     &   phase * sum( zz(1:nv,itp)*cphik(ias:ias+nv-1,it) )
          ! c end of t'(unoccupied)-loop
          !1      end do
          ! c end of t(occupied)-loop
          !1      end do
          !3        call zgemm('T','N',nt0,ntp0,nv,
          !3     &   phase, cphik(ias:ias+nv-1,1:nt0),  nv,
          !3     &          dcmplx(rr(1:nv,1:ntp0),-cc(1:nv,1:ntp0)),  nv,
          !3     &          0d0,
          !3    &   zpsi2b(imdim(iap)-1+i,nctot+1:nctot+nt0,1:ntp0), nt0)
          call zgemm('T','N', nt0,ntp0,nv, &
               phase, cphik(ias,1),nlmto, zz,nlnmx, beta, &
               zppb, nzppb1 )
          ib = imdim(iap)-1+i
          zpsi2b(ib,nctot+1:nctot+nt0,1:ntp0)=zppb(1:nt0,1:ntp0)
          !------------------------------------------------------
       end do !end of optimal product basis-loop
    end do !end of atom-loop
    !      deallocate(rr,cc,iasx)
    deallocate(zz,zppb,iasx)
  end subroutine psi2b_v2

  !------------------------------------------------------------------------------------
  subroutine psicb_v2 (icore,ncore,ntp0,iclass,coskt,sinkt, &
       cphikq, ppb,    nlnmv,nlnmc,mdim, &
       imdim,iatomp, &
       mdimx,nlmto,nbloch,nlnmx,nt,ntp,natom,nclass, &
       nl,nnc, &
       zpsi2b)! rpsi2b,cpsi2b)
    ! written by Ferdi  92.03.17
    ! takao modified at Apr 2002

    ! calculates <psi(k+q,t') | core(k,t) B(R,i)>
    ! for all R
    ! psi(k,t) = S[RLn] b(RLn,k,t)*X(RLn,k)
    !          = S[RLn] b(RLn,k,t)*Phi(RLn,k) + hb(RLn,k,t)*Phidot(RLn,k)
    ! core(k,t)= core states
    ! B(R,i)   = Bloch orthonormal product basis for atom R

    ! <psi(k+q,t') | core(k,t) B(R,i)>
    ! = S[RLn]  b(RLn,k+q,t')^* <Phi(RLn)    |core(k,t) B(R,i)>
    ! + S[RLn] hb(RLn,k+q,t')^* <Phidot(RLn) |core(k,t) B(R,i)>

    ! ntp0       = no. unoccupied states
    ! coskt,sinkt= exp(ik.T)
    ! cphikq  = real and imaginary part of b(k+q).
    !            coefficients of eigenfunctions for argumentationwaves in each MT

    ! icore      = index for core states
    ! ncore      = no. core states in each class
    ! ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>

    ! nlnmv      = number of l,n,m for valence
    ! nlnmc      = number of l,n,m for core states
    ! mdim       = number of optimal product basis functions for each class
    ! nbloch     = total no. optimal product basis
    ! nlnmx      = maximum number of l,n,m
    ! nt         = maximum number of occupied states
    ! ntp        = maximum number of unoccupied states

    ! zpsi2b     =  the matrix elements

    implicit real*8(a-h,o-z)
    implicit integer (i-n)
    complex(8):: cphikq(nlmto,ntp0),zpsi2b(nbloch,nt,ntp),phase
    dimension &
         !                rbkq(nlmto,ntp0),cbkq(nlmto,ntp0),
         !     i          rhbkq(nlmto,ntp0),chbkq(nlmto,ntp0),
         icore(nl*nl*nnc,nclass),ncore(nclass), &
         ppb(nlnmx,nlnmx,mdimx,nclass), &
         !     i          pdb(nlnmx,nlnmx,mdimx,nclass),
         !     i          dpb(nlnmx,nlnmx,mdimx,nclass),
         !     i          ddb(nlnmx,nlnmx,mdimx,nclass),
         nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom), &
         coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)
    !      dimension rpsi2b(nbloch,nt,ntp),
    !     o          cpsi2b(nbloch,nt,ntp)
    ! initialise matrix elements
    !      call dinit   (rpsi2b,nbloch*nt*ntp)
    !      call dinit   (cpsi2b,nbloch*nt*ntp)

    zpsi2b = 0d0
    ! loop over atoms
    ib         = 0
    ias        = 1
    ics        = 0
    do      ia = 1,natom
       ic         = iclass(ia)
       nc         = nlnmc(ic)
       nv         = nlnmv(ic)
       nc1        = nc + 1
       phase  =  dcmplx(coskt(ia),sinkt(ia))
       ! loop over optimal product basis
       iap        = iatomp(ia)
       icp        = iclass(iap)
       ib         = imdim(iap)-1
       do       i = 1,mdim(icp)
          ib         = ib + 1

          ! S[Ln] bkq(Ln,t')^(*) * <Phi(Ln) core(L'n') B(i)> etc.
          ! for a given i, for all L'n' and t'
          ! bkq is complex but < > is real
          do     itp = 1,ntp0
             do      it = 1,ncore(ic)
                icr        = icore(it,ic)

                ! cccccccccccc
                !      write(6,*),it,ic,icore(it,ic)
                ! cccccccccccc

                !      rs1        = vdv ( rbkq(ias,itp),ppb(nc1,icr,i,icp),nv)
                !     .           + vdv (rhbkq(ias,itp),dpb(nc1,icr,i,icp),nv)
                !      cs1        = vdv ( cbkq(ias,itp),ppb(nc1,icr,i,icp),nv)
                !     .           + vdv (chbkq(ias,itp),dpb(nc1,icr,i,icp),nv)

                !      rpsi2b(ib,ics+it,itp) = rs1*coskt(ia) + cs1*sinkt(ia)
                !      cpsi2b(ib,ics+it,itp) = rs1*sinkt(ia) - cs1*coskt(ia)
                ! cccccccccccccccccccccccccccccc
                !      if(abs(rs1)>1d8.or.abs(cs1)>1d8) then
                !        print *,'  psicb2*:'
                !        print *, nc1,icr,i,icp, ppb(nc1,icr,i,icp)
                !        print *,    dpb(nc1,icr,i,icp)
                !       stop
                !      endif
                ! cccccccccccccccccccccccccccccc

                zpsi2b(ib,ics+it,itp) = phase* &
                     dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))

                ! end of t'(unoccupied)-loop
             end do
             ! end of t(occupied)-loop
          end do

          ! end of optimal product basis-loop
       end do

       ! end of atom-loop
       ias        = ias + nlnmv(ic)
       ics        = ics + ncore(ic)
    end do

    return
  end subroutine psicb_v2


  !$$$c -------------------------------------------------------------------
  !$$$      subroutine matzwz2(zw,zmel, ntp0,nstate,ngb, zwz)
  !$$$      implicit none
  !$$$      integer(4) :: nstate,ntp0,itp,it,itp2,it2,ngb
  !$$$      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),
  !$$$     c              zwz(ntp0,nstate,nstate,ntp0)
  !$$$      complex(8), allocatable :: CC(:,:,:)
  !$$$      allocate(CC(ngb,nstate,ntp0) )
  !$$$      call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
  !$$$      do itp2 = 1,ntp0
  !$$$      do  it2 = 1,nstate
  !$$$      do  it  = 1,nstate
  !$$$      do itp  = 1,ntp0
  !$$$        zwz(itp,it,it2,itp2)
  !$$$     & = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
  !$$$      enddo
  !$$$      enddo
  !$$$      enddo
  !$$$      enddo
  !$$$      deallocate(CC)
  !$$$      end
  ! -------------------------------------------------------------------
  subroutine matzwz3(zw,zmel1,zmel2, ntp0,nstate,ngb, zwz)
    implicit none
    integer(4) :: nstate,ntp0,itp,it,itp2,it2,ngb
    complex(8) :: zw(ngb,ngb),zmel1(ngb,nstate,ntp0), &
         zmel2(ngb,nstate,ntp0), &
         zwz(ntp0,nstate,nstate,ntp0)
    complex(8), allocatable :: CC(:,:,:)
    allocate(CC(ngb,nstate,ntp0) )
    call matm(zw,zmel2,cc, ngb, ngb, nstate*ntp0)
    do itp2 = 1,ntp0
       do  it2 = 1,nstate
          do  it  = 1,nstate
             do itp  = 1,ntp0
                zwz(itp,it,it2,itp2) &
                     = sum( dconjg(zmel1(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
             enddo
          enddo
       enddo
    enddo
    deallocate(CC)
  end subroutine matzwz3

end subroutine wmatqk_mpi

logical function newansisoW()
  newansisoW=.true.
end function newansisoW

subroutine readppovl0(q,ngc,ppovl)
  implicit none
  integer, intent(in) :: ngc
  complex(8), intent(out) :: ppovl(ngc,ngc)
  real(8), intent(in) :: q(3)
  integer:: ngc_r,ippovl0
  real(8):: qx(3),tolq=1d-8
  open(newunit=ippovl0,file='PPOVL0',form='unformatted')
  do
     read(ippovl0) qx,ngc_r
     if(sum(abs(qx-q))<tolq) then
        if(ngc_r/=ngc) call rx( 'readin ppovl: ngc_r/=ngc')
        read(ippovl0) ppovl
        exit
     endif
  enddo
  close(ippovl0)
end subroutine readppovl0


