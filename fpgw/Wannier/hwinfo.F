      program hwinfo

c-------------------------------------------------
c Output informations on the cell, atoms, orbitals, ...
c from hpsig
c
c Takashi Miyake, June 2007
c-------------------------------------------------
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


      complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:)
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
     &              ia,iwf,nphix,ib,ifout
      integer(4),allocatable:: ikidx(:),ikbidx(:,:),iphi(:,:),iphidot(:,:),
     &                         nphi(:)
      real(8),allocatable :: bbv(:,:),r0g(:,:),c1(:,:,:),c2(:,:,:),
     &                       phig(:,:,:,:),wphi(:,:)
      real(8) :: pgnorm,wgt,ndg(3),sij,wphis
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
!        write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
      enddo
      endif
      print *,' ======== nqbz qbz  =',nqbz
      print *,' ======== nqibz ngrp=',nqibz,ngrp
!      print *,  qbz 
!      print *,'============================'
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
!      ifhbe      = iopen('hbe.d',1,0,0)
!      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
!      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
!      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'

c --- Readin nlam index
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      ldim2 = nlmto
      read(ifoc,*)
      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
      do ix =1,ldim2
        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
      enddo   

ctm
!      ifbb = iopen('BBVEC',1,0,0)
!      read(ifbb,*)
!      read(ifbb,*)nbb,nqbz2
!      if (nqbz.ne.nqbz2) stop 'readbb: nqbz is wrong!'
!      allocate (bbv(3,nbb),ikbidx(nbb,nqbz))
!      call  readbb(ifbb,nqbz,nspin,nbb,
!     o             bbv,
!     o             ikbidx,
!     o             iko_ixs,iko_fxs,noxs)

c GWinput data
      call getkeyvalue("GWinput","<MLWF>",unit=ifmloc,status=ret) 
      print *,' ifmloc ret=',ifmloc,ret
      read (ifmloc,*) nwf
ctm 040928
      allocate (nphi(nwf))
      read (ifmloc,*) (nphi(i),i=1,nwf)
      nphix = 0
      do i = 1,nwf
         if(nphi(i).gt.nphix)nphix = nphi(i)
      enddo
      allocate (r0g(nphix,nwf),iphi(nphix,nwf),iphidot(nphix,nwf),
     &          wphi(nphix,nwf))
      wphi = 0d0
      do i=1,nwf
         do j=1,nphi(i)
            read(ifmloc,*)iphi(j,i),iphidot(j,i),r0g(j,i),wphi(j,i)
         enddo
      enddo
      close(ifmloc)

c 061004
!      call getkeyvalue("GWinput","wan_gauss_head",ghead,default=.false.)
!      call getkeyvalue("GWinput","wan_truncate",tailt,default=.false.)

c normalize wphi
      do i = 1,nwf
         wphis = dsqrt(sum(wphi(1:nphi(i),i)**2))
         wphi(1:nphi(i),i) = wphi(1:nphi(i),i)/wphis
      enddo

!      write(*,*)'Gaussian Head =',ghead
      write(*,*)'Num. Wannier =',nwf
      write(*,*)'        No.   n    l    m  ibas'
      do i = 1,nwf
         write(*,*)'iwf,',i,nphi(i)
         do j = 1,nphi(i)
            ix = iphi(j,i)
            write(*,"('phi   ',5i5)")
     &      ix,n_indx(ix),l_indx(ix),m_indx(ix),ibas_indx(ix)
            ix = iphidot(j,i)
            write(*,"('phidot',5i5)")
     &      ix,n_indx(ix),l_indx(ix),m_indx(ix),ibas_indx(ix)
         enddo
         write(*,*)
      enddo


      ifout = iopen('wan.info',1,-1,0)
      write(ifout,*)'*** Number of Wannier functions (WFs)'
      write(ifout,*)nwf
      write(ifout,*)'*** Initial guess ***' 
      write(ifout,*)'*** The i-th WF is constructed from N orbitals.'
      do i = 1,nwf
         write(ifout,992)i,nphi(i)
      enddo
      write(ifout,*)'*** The j-th (j=1,2,...,N) component of the i-th WF' 
      write(ifout,*)'   i    j    l    m   ibas    x           y           z          weight'
      do i = 1,nwf
         do j = 1,nphi(i)
            ix = iphi(j,i)
            ib = ibas_indx(ix)
            write(ifout,993)i,j,
     &      l_indx(ix),m_indx(ix),ib,
     &      pos(1,ib),pos(2,ib),pos(3,ib),wphi(j,i)
c            write(*,"('phi   ',5i5)")
c     &      ix,n_indx(ix),l_indx(ix),m_indx(ix),ibas_indx(ix)
c            ix = iphidot(j,i)
c            write(*,"('phidot',5i5)")
c     &      ix,n_indx(ix),l_indx(ix),m_indx(ix),ibas_indx(ix)
         enddo
c         write(*,*)
      enddo
      write(ifout,*)'*** alat'
      write(ifout,*)alat
      write(ifout,*)'*** unit cell' 
      write(ifout,996)plat(1:3,1)
      write(ifout,996)plat(1:3,2)
      write(ifout,996)plat(1:3,3)


 991  format(1x,a39,i5)
 992  format(5x,i5,23x,i5)
 993  format(2i5,3i5,3f12.6,f12.6)
 994  format(1x,a10,f12.6)
 995  format(1x,a33)
 996  format(3f12.6)

      print *, ' ====== end ========================================'
      stop 'hwinfo: ok'
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
      subroutine  readbb(ifbb,nqbz,nspin,nbb,
     o                   bbv,
     o                   ikbidx,
     o                   iko_ixs,iko_fxs,noxs)

      implicit real*8(a-h,o-z)
      parameter (eps = 1d-4)

      real (8) :: u(3),bbv(3,nbb)
      integer (4) :: iopen,
     i           iko_ixs(2),iko_fxs(2),noxs(2)
      integer(4):: ikbidx(nbb,nqbz)
c      integer(4),allocatable:: ikidx(:),ikbidx(:,:)

c      ifbb = iopen('BBVEC',1,0,0)
c      read(ifbb,*)
c      read(ifbb,*)nbb,nqibz2,nqbz2
c      if (nqibz.ne.nqibz2) stop 'readbb: nqibz is wrong!'
c      if (nqbz.ne.nqbz2) stop 'readbb: nqbz is wrong!'
c
c      allocate (ikidx(nqibz),ikbidx(nbb,nqibz))

      do i = 1,nbb
         read(ifbb,*)bbv(1,i),bbv(2,i),bbv(3,i),dummy4
      enddo
      do iq = 1,nqbz
         read(ifbb,*)itmp,u(1:3)
         do ib = 1,nbb
            read(ifbb,*)itmp,itmp2,ikbidx(ib,iq),u(1:3)
         enddo
      enddo
      read(ifbb,*)
      read(ifbb,*)nspin2
      if (nspin .ne. nspin2) stop 'nspin is wrong!'
      do is = 1,nspin
         read(ifbb,*)iko_ixs(is),iko_fxs(is),noxs(is)
      enddo


      return
      end
c-----------------------------------------------------------------------
