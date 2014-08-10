      program hphig

c-------------------------------------------------
c from huumat.f
c Expand |g> by phi; |g> = S[n]c(n)|phi(n)>
c output: coefficients c(n)'s 
c
c g   : c1*phi + c2*phidot (MT) 
c       gaussian           (otherwise)
c phi : atomic orbital (MT)
c       planewave      (otherwise)
c
c Approximation
c in the MTs except for the central MT, phi and phidot are 
c replaced with plane-waves
c
c Takashi Miyake, June 2004
c-------------------------------------------------
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue

      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     & , ngvecpf1(:,:), ngvecpf2(:,:), ngvecp(:,:), 
     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),
     & symope(:,:),phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:)
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
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecb_o,mrecg_o,
     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
     &   ifepscond,nxx,ifvxcpout,ifgb0vec
     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax  !,ifev1,ifev2
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &   ,ngpmx,ngpmx2,ikpisp,ierr

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dini,
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


      complex(8),allocatable:: geig1(:,:),geig2(:,:),
     & cphi1(:,:),cphi2(:,:),uum(:,:,:), ppovl(:,:),
     & ppovl2(:,:)
      complex(8):: ppj,phaseatom
      real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
      integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin
     & ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx
     &, iopen,ico,lxd,lx,ll
      real(8):: ylk

ctm
      integer(4) :: ixc,idummy,idummy2,i1,i2,i3,nbbloop,
     &              ifcphig,ifgeigg,ifdim,
     &              ifpsig(2),ifmloc,ret,nwf,
     &              ifbb,nbb,iko_ixs(2),iko_fxs(2),noxs(2),
     &              iqibz,iqbz,ibb,itmp,itmp2,iti,itf,
     &              nqibz2,nqbz2,iqb,ibb2,iqtmp,ibbtmp,
     &              ia,iwf,nphix,isx,ifiqg,nqtt,nnnn,ikp
      integer(4),allocatable:: ikidx(:),ikbidx(:,:),iphi(:,:),
     &                         iphidot(:,:),nphi(:),ngp(:)
      real(8),allocatable :: bbv(:,:),r0g(:,:),c1(:,:,:),c2(:,:,:),
     &                       phig(:,:),wphi(:,:)
     &                      ,wk(:,:),aau(:,:,:,:)
     &                      ,phiphi(:,:,:),dnorm(:,:)
     &                      ,qtt(:,:)
      real(8) :: pgnorm,wgt,ndg(3),sij,wphis,QpGcut_psi
      complex(8) :: ctmp
      complex(8),allocatable :: psig(:,:,:),cphig(:,:)
     &                         ,qgg(:,:,:),qgg2(:,:,:)
     &                         ,qgg3(:,:),qgg4(:)

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
      call readngmx('QGpsi',ngpmx)
      print *,' ngpmx=',ngpmx

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
      call readngmx('QGpsi',ngpmx)
c      allocate( geig1(ngpmx,nband) )
c      print *, 'end of initialization'

c --- Readin nlam index
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      ldim2 = nlmto
      read(ifoc,*)
      allocate(m_indx(ldim2),n_indx(ldim2),
     &         l_indx(ldim2),ibas_indx(ldim2))
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

c normalize wphi
      do i = 1,nwf
         wphis = dsqrt(sum(wphi(1:nphi(i),i)**2))
         wphi(1:nphi(i),i) = wphi(1:nphi(i),i)/wphis
      enddo

c norm
      allocate(dnorm(nwf,nspin))
      dnorm=0d0

c matching condition at MT radius
      allocate(c1(nphix,nwf,nspin),c2(nphix,nwf,nspin))
      call getc1c2(phitotr,r0g,
     i             iphi,iphidot,n_indx,l_indx,m_indx,ibas_indx,
     i             rr,nrofi,rmax,nphi,
     d             nwf,nlmto,nrx,nl,nn,nbas,nspin,nphix,
     o             c1,c2)

      write(*,*)'# gaussian =',nwf
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

C --- | g_i > = S[a,u]a(a,u;i)*phi(a,u)
c i   : nphix * nwf * nspin (index for local orbital)
c a,u : nlmto (index for atomic orbital phi and phidot)
      allocate ( aau(nlmto,nphix,nwf,nspin)
     &          ,cphig(nlmto,nwf),phig(nlmto,nwf)
     &          ,phiphi(nlmto,nlmto,nspin) )
      call aaumat(iphi,iphidot,c1,c2,phitoto,phitotr,
     i             rr,nrofi,rmax,aa,bb,
     i             n_indx,l_indx,m_indx,ibas_indx,nphi,
     d             nwf,nlmto,nn,nl,nbas,nrx,nspin,nphix,
     o             aau)
c      write(*,*)'aaumat done'
      call ophiphi(phitoto,phitotr,
     i             rr,nrofi,rmax,aa,bb,
     i             n_indx,l_indx,m_indx,ibas_indx,
     d             nwf,nlmto,nn,nl,nbas,nrx,nspin,
     o             phiphi)
c      write(*,*)'ophiphi done'

c file open
      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband
      ifcphig = iopen('CPHIg',0,-1,mrecb_o)
      ifgeigg = iopen('GEIGg',0,-1,mrecg_o)
      ifdim = iopen('PHIG.d',1,-1,0)
      write(ifdim,*) nspin,nqbz,nwf,nlmto,ngpmx
c      ifpsig(1) = iopen('PHIGU',0,-1,0)
c      write(ifpsig(1))nqbz,nwf
c      if (nspin.eq.2) then
c        ifpsig(2) = iopen('PHIGD',0,-1,0)
c        write(ifpsig(2))nqbz,nwf
c      endif

c======================================================================
c --- Set q1(j1range)
c======================================================================
! Note that the true q when we generate eigenfunctions are q1x and q2x.
! q1-q1x should be a G vector.
!
! --- I inserted checkagree to make sure that q1=q1x and q2=q2x ... 
!
c======================================================================

c... Readin QGpsi
      ifiqg  = 101
      open(ifiqg ,file='QGpsi',form='unformatted')
      read(ifiqg) nqtt , ngpmx2, QpGcut_psi, nnnn
      write(6,*)'read(ifiqg)', nqtt , ngpmx, QpGcut_psi, nnnn
      allocate( ngvecp(3,ngpmx),qtt(3,nqtt),ngp(nqtt) )
      do ikp = 1,nqtt
        read (ifiqg) qtt(:,ikp), ngp(ikp) 
        read (ifiqg) ngvecp(1:3, 1:ngp(ikp))
      enddo
      close(ifiqg)
      deallocate(ngvecp)

ctm
c      do 1070 iqbz = 1,nqbz
      do 1070 iqbz = 1,nqtt
         write(*,*)'iqbz =',iqbz, 'out of',nqtt

c         q1(:) = qbz(:,iqbz)
         q1(:) = qtt(:,iqbz)

C --- q1x
      call readqg0('QGpsi',q1,ginv,  q1x, ngp1)
c      call checkagree(q1,q1x,' q1 ne q1x')
      write(6,"('q1 q1x=',3f9.4,3x,3f9.4)") q1,q1x

c ... lxx and allocations
      lxx=2*(nl-1)

C --- < q+G | g >
      allocate( ngvecpf1(3,ngp1),ppovl(ngp1,ngp1)
     &         ,qgg(ngp1,nphix,nwf) )
      call readqg('QGpsi',q1,ginv, q1x, ngp1, ngvecpf1)
c      call checkagree(q1,q1x,' q1 ne q1x xxx2')
      call qggmat(q1,alat,plat,qbas,
     i            iphi,r0g,pos,
     i            rr,nrofi,rmax,aa,bb,
     i            n_indx,l_indx,m_indx,ibas_indx,
     i            ngvecpf1,nphi,
     d            lxx,ngp1,nwf,nlmto,nbas,nrx,nphix,
     o            qgg)
c      write(*,*)'qggmat done'
      call mkppovl2(alat,plat,qbas, !22April2004 
     &    ngp1, ngvecpf1, 
     &    ngp1, ngvecpf1, 
     &    nbas, rmax, pos,  
     o    ppovl)

c debug
c      allocate (qgg4(ngp1))
c
      allocate (qgg2(ngp1,nphix,nwf),ppovl2(ngp1,ngp1))
      ppovl2=ppovl
      qgg2=0d0
      call matcinv(ngp1,ppovl2)
      do iwf = 1,nwf
      do j   = 1,nphi(iwf)
         call matm(ppovl2,qgg(:,j,iwf),qgg2(:,j,iwf),ngp1,ngp1,1)
c debug
c         call matm(ppovl,qgg2(:,j,iwf),qgg4(:),ngp1,ngp1,1)
c         qgg4 = qgg4 - qgg(:,j,iwf)
c         call matm(ppovl2,qgg(:,j,iwf),qgg4(:),ngp1,ngp1,1)
c         qgg4 = qgg4 - qgg2(:,j,iwf)
c         do i=1,ngp1
c           if ((dabs(real(qgg4(i))).gt.1d-8).or.
c     &         (dabs(imag(qgg4(i))).gt.1d-8)) then
c           write(*,*)iwf,j,i
c           write(*,*)qgg4(i)
c           endif
c         enddo
      enddo
      enddo
      deallocate (ppovl2)
c      deallocate (qgg4)

C --- Calcuate <psi{q1x j1} | g_{iwf}> 
      do 1050 ispin=1,nspin

      if (iqbz.eq.1) then
c ... MT part
cr   ldim2 = nlmto
cr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
cr   l_indx   (1:ldim2) : l index
cr   ibas_indx(1:ldim2) : ibas index.
       phig = 0d0
       do 1020 ia = 1,nlmto
       do iwf = 1,nwf
       do j   = 1,nphi(iwf)
         phig(ia,iwf) =  phig(ia,iwf)
     &      + wphi(j,iwf)*aau(ia,j,iwf,ispin)
       enddo
       enddo
 1020  continue
c       write(*,*)'MT part done'

       do iwf=1,nwf
       do ia2=1,nlmto
       do ia1=1,nlmto
          dnorm(iwf,ispin) = dnorm(iwf,ispin) 
     &    + phig(ia1,iwf)*phig(ia2,iwf)*phiphi(ia1,ia2,ispin)
       enddo
       enddo
       enddo

       cphig(:,:)=phig(:,:)
       deallocate(phig)
c end of if iqbz=1
      endif

c ... Interstitial Plane Wave part
       allocate( qgg3(ngpmx,nwf) )
c       call readgeig(q1, ngpmx, ispin, quu, geig1)
c       call checkagree(q1,quu,' q1 ne quu eig')
       qgg3 = 0d0
       do iwf = 1,nwf
       do j   = 1,nphi(iwf)
          qgg3(1:ngp1,iwf) = qgg3(1:ngp1,iwf)
     &  +  wphi(j,iwf)* qgg2(1:ngp1,j,iwf)
       enddo
       enddo

      if (iqbz.le.nqbz) then
       do iwf=1,nwf
       do ig2=1,ngp1
       do ig1=1,ngp1
          dnorm(iwf,ispin) = dnorm(iwf,ispin) 
     &    + dcmplx(qgg3(ig1,iwf))*qgg3(ig2,iwf)
     &    *ppovl(ig1,ig2) / dble(nqbz)
       enddo
       enddo
       enddo
      endif

c       deallocate(cphi1)


ctm
         ikpisp = ispin + nspin*(iqbz-1)
         write(ifgeigg, rec=ikpisp) qgg3(1:ngpmx,1:nwf)
         write(ifcphig, rec=ikpisp) cphig(1:ldim2,1:nwf)

c         do iwf=1,nwf
c         do i=1,ldim2
c            write(*,*)i,cphig(i,iwf)
c         enddo
c         enddo
 1050 continue

      deallocate(ngvecpf1, ppovl, qgg, qgg2, qgg3)

c end of qbz-loop
 1070 continue


c dnorm = |<phi|phi>|^2 --> dnorm = |<phi|phi>|^(1/2)
      do ispin=1,nspin
      do iwf=1,nwf
         dnorm(iwf,ispin)=dnorm(iwf,ispin)**(0.25d0)
         write(*,*)'norm: ',ispin,iwf,dnorm(iwf,ispin)**2
      enddo
      enddo

c renormalize
      allocate(qgg3(ngpmx,nwf))
c      do iqbz=1,nqbz
      do iqbz=1,nqtt
      do ispin=1,nspin
         ikpisp = ispin + nspin*(iqbz-1)
         read(ifgeigg, rec=ikpisp) qgg3(1:ngpmx,1:nwf)
         read(ifcphig, rec=ikpisp) cphig(1:ldim2,1:nwf)
         do iwf=1,nwf
            qgg3(:,iwf)  = qgg3(:,iwf)/dnorm(iwf,ispin)
            cphig(:,iwf) = cphig(:,iwf)/dnorm(iwf,ispin)
         enddo
         write(ifgeigg, rec=ikpisp) qgg3(1:ngpmx,1:nwf)
         write(ifcphig, rec=ikpisp) cphig(1:ldim2,1:nwf)
      enddo
      enddo
      deallocate(qgg3,cphig)

      print *, ' ====== end ========================================'
c      stop ' ====== end ========================================'
      stop 'hpsig: ok'
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
      subroutine  getc1c2(phitotr,r0g,
     i             iphi,iphidot,n_indx,l_indx,m_indx,ibas_indx,
     i             rr,nrofi,rmax,nphi,
     d             nwf,nlmto,nrx,nl,nn,nbas,nspin,nphix,
     o             c1,c2)

      implicit real*8(a-h,o-z)
      real(8) :: phitotr(nrx,0:nl-1,nn,nbas,nspin),r0g(nphix,nwf),
     .           rr(nrx,nbas),rmax(nbas),
     .           c1(nphix,nwf,nspin),c2(nphix,nwf,nspin)
      integer(4) :: iphi(nphix,nwf),iphidot(nphix,nwf),
     .              n_indx(nlmto),l_indx(nlmto),m_indx(nlmto),
     .              ibas_indx(nlmto),nrofi(nbas),nphi(nwf)


      pi = 4d0*datan(1d0)

      c1 = 0d0
      c2 = 0d0

      do is  = 1,nspin
      do iwf = 1,nwf
      do j   = 1,nphi(iwf)

         ip = iphi(j,iwf)
         ibp = ibas_indx(ip)
         inp = n_indx(ip)
         ilp = l_indx(ip)
         imp = m_indx(ip)

         id = iphidot(j,iwf)
         ibd = ibas_indx(id)
         ind = n_indx(id)
         ild = l_indx(id)
         imd = m_indx(id)

         if (ibp .ne. ibd) stop 'getc1c2: ibp.ne.ibd'
         if (ilp .ne. ild) stop 'getc1c2: ilp.ne.ild'
         if (imp .ne. imd) stop 'getc1c2: imp.ne.imd'

c g(r) = A exp(-(r/r0)**2)
c A = (128 / (pi*r0))**(1/4)
         Ag = (128d0 / (pi*r0g(j,iwf)**6))**0.25d0
         rmt = rr(nrofi(ibp),ibp)
         g  = Ag * dexp(-(rmt/r0g(j,iwf))**2)
         dg = g * (-2d0*rmt/r0g(j,iwf)**2)

c phi
         rmt  = rr(nrofi(ibp),ibp)
         rmt1 = rr(nrofi(ibp)-1,ibp)
         p  =  phitotr(nrofi(ibp),ilp,inp,ibp,is)/rmt
         dp = (phitotr(nrofi(ibp),  ilp,inp,ibp,is)/rmt
     &       - phitotr(nrofi(ibp)-1,ilp,inp,ibp,is)/rmt1)
     &       / (rmt - rmt1)

c phidot
         d  =  phitotr(nrofi(ibd),ild,ind,ibd,is)/rmt
         dd = (phitotr(nrofi(ibd),  ild,ind,ibd,is)/rmt
     &       - phitotr(nrofi(ibd)-1,ild,ind,ibd,is)/rmt1)
     &       / (rmt - rmt1)

c c1 *  p + c2 *  d =  g
c c1 * dp + c2 * dd = dg
         detinv = 1d0 / (p * dd - dp * d)
         c1(j,iwf,is) = detinv * (dd  * g - d * dg)
         c2(j,iwf,is) = detinv * (-dp * g + p * dg)

check
          gtmp = c1(j,iwf,is) *  p + c2(j,iwf,is) *  d
         dgtmp = c1(j,iwf,is) * dp + c2(j,iwf,is) * dd
         write(*,*)'iwf,j =',iwf,j
         write(*,*)'c1,c2             ',c1(j,iwf,is),c2(j,iwf,is)
         write(*,*)'phi,phidot        ',p,d
         write(*,*)'dphi/dr,dphidot/dr',dp,dd
         write(*,*)'g,gtmp            ',g,gtmp
         write(*,*)'dg,dgtmp          ',dg,dgtmp

      enddo
      enddo
      enddo

      return
      end
c--------------------------------------------------------------
      subroutine qggmat(q,alat,plat,qlat,
     i            iphi,r0g,bas,
     i            rr,nrofi,rmax,aa,bb,
     i            n_indx,l_indx,m_indx,ibas_indx,
     i            ngvec,nphi,
     d            lxx,ng,nwf,nlmto,nbas,nrx,nphix,
     o            qgg)
c < q+G | g> where q+G denotes IPW, zero within MT sphere
c g(r) = Ag exp(-(r/r0)**2)
      implicit real*8 (a-h,o-z)

      parameter (nmt = 4)
      parameter (eps = 1d-5)

      complex(8) :: qgg(ng,nphix,nwf),eiqgr
      complex(8) :: img =(0d0,1d0)
      real(8) :: q(3),plat(3,3),qlat(3,3),r0g(nphix,nwf),bas(3,nbas),
     &  rr(nrx,nbas),rmax(nbas),aa(nbas),bb(nbas),
     &  g(3),qg(3),ag(nphix,nwf)
      integer(4) :: iphi(nphix,nwf),nrofi(nbas),ntail(nphix,nwf),
     &  n_indx(nlmto),l_indx(nlmto),m_indx(nlmto),ibas_indx(nlmto),
     &  ngvec(3,ng),nphi(nwf)

      integer(4) :: ng,nwf,nlmto,nbas,nrx
      real(8) :: absqg,tripl,pi,rbas(3,nphix,nwf),
     &  alat,tpiba,tpibaqlat(3,3),voltot 
      real(8),allocatable :: cy(:),yl(:),rtail(:,:,:),rg(:,:,:),
     &                       phij(:),psij(:),rjl(:)

c-----------------------------------------------------
Cr True q is given by
Cr    True_q(1:3)     = 2*pi/alat * q(1:3)
Cr True G is given by 
Cr    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvec(1:3,ig)) ,ig=1,ng
c-----------------------------------------------------

c      print *,' qggmat:'
      pi        = 4d0*datan(1d0)
      voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      tpibaqlat =  2*pi/alat *qlat
      tpiba     =  2*pi/alat

      qgg = 0d0
      ntailx = 0
      do iwf = 1,nwf
      do j   = 1,nphi(iwf)
         ilmto = iphi(j,iwf)
         ibas  = ibas_indx(ilmto)
         rbas(1:3,j,iwf) =  bas(1:3,ibas)*alat
         ag(j,iwf) = (128d0 / (pi*r0g(j,iwf)**6))**0.25d0
         ntail(j,iwf) = int(log(dble(nmt)) / aa(ibas))
     &              + nrofi(ibas) + 1
         if (ntailx .lt. ntail(j,iwf)) ntailx = ntail(j,iwf)
      enddo
      enddo

      allocate(cy((lxx+1)**2),yl((lxx+1)**2),
     &         rtail(ntailx,nphix,nwf),rg(ntailx,nphix,nwf),
     &         phij(0:lxx),psij(0:lxx),rjl(ntailx))


      do iwf = 1,nwf
      do j   = 1,nphi(iwf)
         ilmto = iphi(j,iwf)
         ibas =  ibas_indx(ilmto)
         do ir = 1,ntail(j,iwf)
            rtmp = bb(ibas)*exp(aa(ibas)*dble(ir-1))
            rtail(ir,j,iwf) = rtmp
c rg = r * gaussian
            rg(ir,j,iwf) = rtmp*ag(j,iwf)*dexp(-(rtmp/r0g(j,iwf))**2)
         enddo
      enddo
      enddo


      do ig  = 1, ng
        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvec(1:3,ig)))
        absqg2  = sum(qg**2)
        absqg   = sqrt(absqg2)

      if (absqg .gt. eps) then

        do iwf = 1,nwf
        do j   = 1,nphi(iwf)

           ilmto = iphi(j,iwf)
           ibas  = ibas_indx(ilmto)
           il    = l_indx(ilmto)
           im    = m_indx(ilmto)
           eiqgr = exp( -img* sum(qg(1:3)*rbas(1:3,j,iwf)) )
           a     = aa(ibas)
           b     = bb(ibas)

c Ylm(qg)
           call sylmnc(cy,lxx)
           call sylm(qg/absqg,yl,lxx,r2s) !spherical factor Y(dq)
           ilm   = il**2 + il+1 + im
           ylmqg = cy(ilm)*yl(ilm)

c r \times j_l(|qg|r)  !bessel function
           do ir =1,ntail(j,iwf)
             rtmp = rtail(ir,j,iwf)
             call bessl(absqg2*rtmp**2,lxx,phij,psij)
c  phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
             if(rtmp==0d0) then
               rjl(ir)=0d0
             else   
               rjl(ir) = rtmp* phij(il)* (absqg*rtmp)**il  
             endif  
           enddo

c I[rmt:\infty] dr r^2 jl(|qg|r) gaussian(r)
           nr = nrofi(ibas)
           call gintxx(rg(1:nr,j,iwf),rjl(1:nr),a,b,nr,sum1)
           nr = ntail(j,iwf)
           call gintxx(rg(1:nr,j,iwf),rjl(1:nr),a,b,nr,sum2)
           r2jlg = sum2 - sum1

c < qg | gaussian >
           qgg(ig,j,iwf) = 4d0*pi*(-img)**il * ylmqg * r2jlg*eiqgr

        enddo
c end of iwf-loop
        enddo

c if (absqg=0)
      else

        do iwf = 1,nwf
        do j   = 1,nphi(iwf)
           ilmto = iphi(j,iwf)
           ibas  = ibas_indx(ilmto)
           il    = l_indx(ilmto)
           im    = m_indx(ilmto)
           a     = aa(ibas)
           b     = bb(ibas)

           if (il.ne.0) then
              qgg(ig,j,iwf) = 0d0
           else
              nr = nrofi(ibas)
              call gintxx(rg(1:nr,j,iwf),rtail(1:nr,j,iwf),a,b,nr,sum1)
              nr = ntail(j,iwf)
              call gintxx(rg(1:nr,j,iwf),rtail(1:nr,j,iwf),a,b,nr,sum2)
              qgg(ig,j,iwf) = (sum2-sum1) * dsqrt(4d0*pi)
           endif
        enddo
        enddo

c end of if (absqg=0)
      endif

c end of ig-loop
      enddo

      deallocate(cy,yl,rtail,rg,phij,psij,rjl)

c      print *,' qggmat: done '
      end

c--------------------------------------------------------------
      subroutine aaumat(iphi,iphidot,c1,c2,phitoto,phitotr,
     i             rr,nrofi,rmax,aa,bb,
     i             n_indx,l_indx,m_indx,ibas_indx,nphi,
     d             nwf,nlmto,nn,nl,nbas,nrx,nspin,nphix,
     o             aau)
c |g> = S[a,u]a(a,u)*phi(a,u)
c g = c1*phi + c2*phidot, which is conneted to the gaussian tail
c     at r=rmt

      implicit real*8 (a-h,o-z)

      real(8) :: aau(nlmto,nphix,nwf,nspin),
     &           phitoto(nrx,0:nl-1,nn,nbas,nspin),
     &           phitotr(nrx,0:nl-1,nn,nbas,nspin)
      real(8) :: rr(nrx,nbas),rmax(nbas),aa(nbas),bb(nbas),
     &           c1(nphix,nwf,nspin),c2(nphix,nwf,nspin),
     &           rphi(nrx),rhead(nrx)
      integer(4) :: iphi(nphix,nwf),iphidot(nphix,nwf),nrofi(nbas),
     &  n_indx(nlmto),l_indx(nlmto),m_indx(nlmto),ibas_indx(nlmto),
     &  nphi(nwf)


      aau = 0d0

      do is  = 1,nspin
      do iwf = 1,nwf
      do j   = 1,nphi(iwf)
         
c r * phi
         inlmp = iphi(j,iwf)
         ibasp = ibas_indx(inlmp)
         inp   = n_indx(inlmp)
         ilp   = l_indx(inlmp)
         imp   = m_indx(inlmp)

c r * phidot
         inlmd = iphidot(j,iwf)
         ibasd = ibas_indx(inlmd)
         ind    = n_indx(inlmd)
         ild    = l_indx(inlmd)
         imd    = m_indx(inlmd)

         if (ibasp .ne. ibasd) stop 'aaumat: bas(phi) .ne. bas(phidot)'
         if (ilp .ne. ild) stop 'aaumat: l(phi) .ne. l(phidot)'
         if (imp .ne. imd) stop 'aaumat: m(phi) .ne. m(phidot)'

         aau(inlmp,j,iwf,is) = c1(j,iwf,is)
         aau(inlmd,j,iwf,is) = c2(j,iwf,is)

      enddo
c end of iwf-loop
      enddo
c end of is-loop
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ophiphi(phitoto,phitotr,
     i             rr,nrofi,rmax,aa,bb,
     i             n_indx,l_indx,m_indx,ibas_indx,
     d             nwf,nlmto,nn,nl,nbas,nrx,nspin,
     o             phiphi)
c < phi | phi >

      implicit real*8 (a-h,o-z)

      real(8) :: phiphi(nlmto,nlmto,nspin),
     &           phitoto(nrx,0:nl-1,nn,nbas,nspin),
     &           phitotr(nrx,0:nl-1,nn,nbas,nspin)
      real(8) :: rr(nrx,nbas),rmax(nbas),aa(nbas),bb(nbas),
     &           rphi(nrx),rhead(nrx)
      integer(4) :: nrofi(nbas),
     &  n_indx(nlmto),l_indx(nlmto),m_indx(nlmto),ibas_indx(nlmto)

      phiphi = 0d0

      do is  = 1,nspin
         
      do inlmp=1,nlmto
         ibasp = ibas_indx(inlmp)
         inp   = n_indx(inlmp)
         ilp   = l_indx(inlmp)
         imp   = m_indx(inlmp)

      do inlmd=1,nlmto
         ibasd = ibas_indx(inlmd)
         ind    = n_indx(inlmd)
         ild    = l_indx(inlmd)
         imd    = m_indx(inlmd)

         if (ibasp .ne. ibasd) cycle
         if (ilp .ne. ild) cycle
         if (imp .ne. imd) cycle

         call gintxx(phitoto(1:nrofi(ibasp),ilp,inp,ibasp,is),
     &               phitoto(1:nrofi(ibasd),ild,ind,ibasd,is),
     &               aa(ibasp),bb(ibasp),nrofi(ibasp),sum)
         phiphi(inlmp,inlmd,is) = sum
c         write(*,*)inlm,iwf,sum
      enddo

      enddo   

c end of is-loop
      enddo

      return
      end
c-----------------------------------------------------------------------
