      program hwfchrac

c-------------------------------------------------
c from hpsig.F
c Calculate <chi|psi>
c
c Takashi Miyake, Mar 2008
c-------------------------------------------------
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue

      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     & , ngvecpf1(:,:), ngvecpf2(:,:), 
     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:),
     &phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:)
      complex(8),parameter:: img=(0d0,1d0)
c,nw,incwf,natom,nclass,ipos,igrp,
c     & iinvg,nspin,nl,nn,nnv,nnc,
c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions 
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
c     o                   izdummy, 
c     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
c     o   iecore,ikonf,iicore,incore,nctot,             !core
c     o   imagw_dummy,niw,idummy,
      integer(4)
     &   nw_input,
     &   ifhbe,
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
     &   ifepscond,nxx,ifvxcpout,ifgb0vec
     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax  !,ifev1,ifev2
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &   ,ngpmx

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini,
     & dwry,dwh,omg_c,omg2,xxx
      integer(4)::nwin, incwfin,  verbose
      real(8)::efin
      integer(4):: bzcase, mrecg,ifphi,
     & nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc,
     & ldim2,ixx,ngp1,ngp2,nq0it
      real(8):: qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
      real(8),parameter::  pi =     3.14159265358979323846d0
      real(8),parameter::  fpi =    4d0*pi

c      logical:: test_qbzreg
      logical:: qbzreg
!-------------------------------------------------------------------------
      integer(4),allocatable:: ncindx(:,:),
     &           lcindx(:,:),
     &           nrad(:),
     &           nindx_r(:,:),
     &           lindx_r(:,:),
     &           nc_max(:,:),
     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
      real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:)
     &                     ,phitotr(:,:,:,:,:),
     &        bb(:),zz(:),rmax(:),cy(:),yl(:)


      complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi(:,:)
     & ,uum(:,:,:), ppovl(:,:)
      complex(8):: ppj,phaseatom
      real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
      integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin
     & ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx
     &, iopen,ico,lxd,lx,ll
      real(8):: ylk

ctm
      integer(4) :: ixc,idummy,idummy2,i1,i2,i3,nbbloop,
     &              ifpsig(2),ifmloc,ret,nwf,
     &              ifbb,nbb,iko_ixs(2),iko_fxs(2),noxs(2),
     &              iqibz,iqbz,ibb,itmp,itmp2,iti,itf,
     &              nqibz2,nqbz2,iqb,ibb2,iqtmp,ibbtmp,
     &              ia,iwf,nphix,iband,isig,nmax
      integer(4),allocatable:: ikidx(:),ikbidx(:,:),iphi(:,:),iphidot(:,:),
     &                         nphi(:),isort(:)
      real(8),allocatable :: bbv(:,:),r0g(:,:),c1(:,:,:),c2(:,:,:),
     &                       phig(:,:,:,:),wphi(:,:)
      real(8) :: pgnorm,wgt,ndg(3),sij,wphis,swphi
      complex(8),allocatable :: psig(:,:,:),qgg(:,:,:)
      logical :: ghead,tailt
c-------------------------

ctm
c mode switch. --------------
c      write(6,*) ' --- Choose modes below ----------------'
c      write(6,*) '  ????????????????????? ' 
c      write(6,*) ' --- Put number above ! -----------------'
c      call readin5(ixc,idummy,idummy2)
c      write(6,*) ' ixc=',ixc
c      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'


c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz  
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()

c--- Use regular mesh even for bzcase==2
      if(bzcase()==2.and.qbzreg()) then
      deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
      do i=1,nqbz
        qbz(:,i) = qbz(:,i) -deltaq/2d0
        write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
      enddo
      endif
      print *,' ======== nqbz qbz  =',nqbz
      print *,' ======== nqibz ngrp=',nqibz,ngrp
      print *,  qbz 
      print *,'============================'
      print *

C--- readin GWIN and LMTO, then allocate and set datas.
      nwin = -999    !readin condition. Not readin NW file
      incwfin= 0     !readin condition. use ForX0 for core in GWIN
      efin =  -999d0 !readin condition. Not readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c---  These are allocated and setted by genallcf_v3
c      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr   
c      character(120):: symgrp
c      character(6),allocatable :: clabl(:)
c      integer(4),allocatable:: iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
c     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c     &    occv(:,:,:),unoccv(:,:,:)   
c     &   ,occc(:,:,:),unoccc(:,:,:),
c     o    nocc(:,:,:),nunocc(:,:,:)
c      real(8), allocatable:: 
c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)

!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) stop ' nclass /= natom '

c --- read dimensions of h,hb
      ifhbe      = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'

c --- read by rdpp ; Radial integrals ppbrd and plane wave part
      call getsrdpp2(nclass,nl,nxx)

c --- read radial functions PHIVC   (taken from hasfp0)
      print *,' Go to readining phivc'
      ifphi  = iopen('PHIVC', 0,-1,0)     ! PHIV+PHIC augmentation wave and core
      read(ifphi) nbas, nradmx, ncoremx,nrx
      if( nbas/=natom ) stop ' nbas(PHIVC) /= natom '
      deallocate(ncore)
      allocate(  ncindx(ncoremx,nbas),
     &           lcindx(ncoremx,nbas),
     &           nrad(nbas),
     &           nindx_r(1:nradmx,1:nbas),
     &           lindx_r(1:nradmx,1:nbas),
     &        aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) ,
     &        phitoto(nrx,0:nl-1,nn,nbas,nspin),
     &        phitotr(nrx,0:nl-1,nn,nbas,nspin),
     &        nc_max(0:nl-1,nbas),ncore(nbas),rmax(nbas) )
      print *,' end of allocation'
      read(ifphi) nrad(1:nbas)
      read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
      nc_max=0
      do ibas=1,nbas
        ic = ibas
        write(6,*)' --- read PHIVC of ibas nrad=',ibas,nrad(ic)
        read(ifphi) ncore(ic), ncoremx                            !core
        read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
        print *,' xxx0'
        read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)

        print *, 'ic icx=',ic,icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
        if(ic/=icx) then
           stop ' h_uu: ic/=icx'
        endif
        print *,' xxx1 ncoremx ncore(ic)=',ncoremx,ncore(ic)
        read(ifphi) rr(1:nrofi(ic),ic)
        print *,' xxx2 ncoremx ncore(ic)=',ncoremx,ncore(ic)

        print *,' xxx2 nspin=',nspin
        rmax(ic) = rr(nrofi(ic),ic)
        do isp = 1, nspin
          write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
          do ico = 1, ncore(ic) !core
            l =  lcindx(ico,ic)
            n =  ncindx(ico,ic)
            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
            phitotr(1:nrofi(ic),l,n, ic,isp)=              !core raw= core orthgonal
     &      phitoto(1:nrofi(ic),l,n, ic,isp)               !
            if(n>nc_max(l,ic)) nc_max(l,ic)=n
            print *,' sss1c=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
          enddo
          do irad = 1, nrad(ic)   !valence
            l = lindx_r (irad,ic)
            n = nindx_r (irad,ic) + nc_max(l,ic)
            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal 
            read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
            print *,' sss1=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
            print *,' sss2=',sum(abs(phitotr(1:nrofi(ic),l,n, ic,isp)))
          enddo
        enddo  
      enddo   

c--- cg coefficient.  y = cg y y ; y is the real spherical harmonics 
      ngrpx=1
      allocate( cg(nl**2,nl**2,(2*nl-1)**2), symope(3,3) )
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
      cg = 0d0 !for sanity check
      call rotcg(nl-1,symope,ngrpx,cg)

c --- initiallization to get eigenfunctions
      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
      call init_readeigen2(mrecb,nlmto,mrecg)

c --- Readin nlam index
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      ldim2 = nlmto
      read(ifoc,*)
      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
      do ix =1,ldim2
        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
      enddo   

c ---  q near zero
      print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
!      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
      nq0it = nq0i
      allocate( q0i(1:3,1:nq0i) ) !wqt(1:nq0i),
!      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      nq0ix = nq0i
      do i=1,nq0i
      read (101,* ) xxx,q0i(1:3,i)
      if(xxx==0d0 ) nq0ix = i-1
      enddo
      nq0i = nq0ix ! New nq0i July 2001
      write(6,*) ' Used k number in Q0P =', nq0i
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
      close(101)


c======================================================================
c --- Set iq,ib,is
      write(*,*)'ispin=',nspin
      read(*,*)ispin
      write(*,*)'iqibz=',nqibz
      read(*,*)iqibz
      write(*,*)'iband=',nband
      read(*,*)iband
c======================================================================

      open(unit=99,file='WFch')
!      do 1070 iqibz = 1,nqibz
!         write(*,*)'iqibz =',iqibz, 'out of',nqibz

         q1(:) = qibz(:,iqibz)
c ... MT part
cr   ldim2 = nlmto
cr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
cr   l_indx   (1:ldim2) : l index
cr   ibas_indx(1:ldim2) : ibas index.
       allocate( cphi (nlmto,nband),wphi(nlmto,nband) )
       call readcphi(q1, nlmto, ispin, quu, cphi)
       wphi = abs(cphi)
c       call checkagree(q1,q1x,' q1 ne quu')

       write(99,*)'***** '
       write(99,"('is=',i5)")ispin
       write(99,"('iq=',i5,', q=',3f12.6)")iqibz,q1
       write(99,"('ib=',i5)")iband
       write(99,*)'***  '
       write(99,*)
     & '   i    m    n    l   bas   abs(cphi)   Re(cphi)    Im(cphi)'
      swphi = sum(wphi(:,iband))
!      cphi = cphi / swphi
!      wphi = wphi / swphi

      do 1020 ia = 1,nlmto
          write(99,"(5i5,3f12.6)")ia,m_indx(ia),n_indx(ia),l_indx(ia),
     &              ibas_indx(ia),wphi(ia,iband),cphi(ia,iband)
 1020  continue

c sorting
       allocate(isort(nlmto))
       call sortea(wphi(:,iband),isort,nlmto,isig)

       write(99,*)'***** Large components'
       write(99,*)
     & '   i    m    n    l   bas   abs(cphi)   Re(cphi)    Im(cphi)'
      nmax = min(12,nlmto)
      do 1030 ia = 1,nmax
          ia2 = isort(nlmto+1-ia)
          write(99,"(5i5,3f12.6)")ia2,m_indx(ia2),n_indx(ia2),
     &       l_indx(ia2),ibas_indx(ia2),wphi(ia2,iband),cphi(ia2,iband)
 1030  continue

       deallocate(cphi,wphi,isort)
c end of qbz-loop
 1070 continue

      stop 'hwfch: ok'
      end
c====================================================================
      subroutine checkagree(a,b,char)
      real(8):: a(3),b(3) 
      character*(*) :: char
      if(sum(abs(a-b))>1d-6) then
         print *,' Error in checkagree:',char
         stop ' Error in checkagree:'
      endif
      end
c-----------------------------------------------------------------------
