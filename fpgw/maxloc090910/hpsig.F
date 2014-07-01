      program hpsig

c-------------------------------------------------
c from huumat.f
c Calculate <psi|g>
c psi : Kohn Sham orbital
c g   : c1*phi + c2*phidot (MT) 
c       gaussian           (otherwise)
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
     &              ia,iwf,nphix
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
      allocate( geig1(ngpmx,nband) )
      print *, 'end of initialization'

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

ctm
      ifbb = iopen('BBVEC',1,0,0)
      read(ifbb,*)
      read(ifbb,*)nbb,nqbz2
      if (nqbz.ne.nqbz2) stop 'readbb: nqbz is wrong!'
      allocate (bbv(3,nbb),ikbidx(nbb,nqbz))
      call  readbb(ifbb,nqbz,nspin,nbb,
     o             bbv,
     o             ikbidx,
     o             iko_ixs,iko_fxs,noxs)

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
      call getkeyvalue("GWinput","wan_gauss_head",ghead,default=.false.)
      call getkeyvalue("GWinput","wan_truncate",tailt,default=.false.)

c normalize wphi
      do i = 1,nwf
         wphis = dsqrt(sum(wphi(1:nphi(i),i)**2))
         wphi(1:nphi(i),i) = wphi(1:nphi(i),i)/wphis
      enddo

c matching condition at MT radius
      allocate(c1(nphix,nwf,nspin),c2(nphix,nwf,nspin))
      if (.not.ghead)
     >call getc1c2(phitotr,r0g,
     i             iphi,iphidot,n_indx,l_indx,m_indx,ibas_indx,
     i             rr,nrofi,rmax,nphi,nc_max,
     d             nwf,nlmto,nrx,nl,nn,nbas,nspin,nphix,
     o             c1,c2)

      write(*,*)'Gaussian Head =',ghead
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

c file open
      ifpsig(1) = iopen('PSIGU',0,-1,0)
      write(ifpsig(1))nqbz,iko_ixs(1),iko_fxs(1),nwf
      if (nspin.eq.2) then
        ifpsig(2) = iopen('PSIGD',0,-1,0)
        write(ifpsig(2))nqbz,iko_ixs(2),iko_fxs(2),nwf
      endif

c======================================================================
c --- Set q1(j1range)
c======================================================================
! Note that the true q when we generate eigenfunctions are q1x and q2x.
! q1-q1x should be a G vector.
!
! --- I inserted checkagree to make sure that q1=q1x and q2=q2x ... 
!
ctm
      j1min = iko_ixs(1)
      j1max = iko_fxs(1)
      if (nspin .eq. 2) then
         if (iko_ixs(2) .lt. j1min) j1min = iko_ixs(2)
         if (iko_fxs(2) .gt. j1max) j1max = iko_fxs(2)
      endif
      j2min = 1
      j2max = nwf
c======================================================================

      allocate( psig(j1min:j1max,nwf,nspin) )

ctm
      do 1070 iqbz = 1,nqbz
         write(*,*)'iqbz =',iqbz, 'out of',nqbz

         q1(:) = qbz(:,iqbz)

C --- q1x
      call readqg0('QGpsi',q1,ginv,  q1x, ngp1)
c      call checkagree(q1,q1x,' q1 ne q1x')
      write(6,"('q1 q1x=',3f9.4,3x,3f9.4)") q1,q1x

c ... lxx and allocations
      lxx=2*(nl-1)

C --- < q+G | g >
      allocate( ngvecpf1(3,ngp1),qgg(ngp1,nphix,nwf) )
      if (tailt) then
        qgg = 0d0
      else
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
      endif

C --- < phi | g >
      allocate ( phig(nlmto,nphix,nwf,nspin) )
      call phigmat(iphi,iphidot,c1,c2,phitoto,phitotr,
     i             rr,nrofi,rmax,aa,bb,r0g,ghead,
     i             n_indx,l_indx,m_indx,ibas_indx,nphi,nc_max,
     d             nwf,nlmto,nn,nl,nbas,nrx,nspin,nphix,
     o             phig)
c      write(*,*)'phigmat done'

C --- Calcuate <psi{q1x j1} | g_{iwf}> 
      psig = 0d0
      do 1050 ispin=1,nspin

c ... MT part
cr   ldim2 = nlmto
cr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
cr   l_indx   (1:ldim2) : l index
cr   ibas_indx(1:ldim2) : ibas index.
       allocate( cphi1 (nlmto,nband) )
       call readcphi(q1, nlmto, ispin, quu, cphi1)
c       call checkagree(q1,q1x,' q1 ne quu')

       do 1020 ia = 1,nlmto
c          do j1= j1min,j1max
          do j1= iko_ixs(ispin),iko_fxs(ispin)
          do iwf = 1,nwf
          do j   = 1,nphi(iwf)
            psig(j1,iwf,ispin) =  psig(j1,iwf,ispin) 
     &       + wphi(j,iwf)*dconjg(cphi1(ia,j1))*phig(ia,j,iwf,ispin)
          enddo
          enddo
          enddo
 1020  continue

c       write(*,*)'MT part done'

c ... Interstitial Plane Wave part
       call readgeig(q1, ngpmx, ispin, quu, geig1)
c       call checkagree(q1,quu,' q1 ne quu eig')
c       do j1=j1min,j1max
       do j1  = iko_ixs(ispin),iko_fxs(ispin)
       do iwf = 1,nwf
       do j   = 1,nphi(iwf)
          psig(j1,iwf,ispin) =  psig(j1,iwf,ispin) 
     &  +  wphi(j,iwf)*sum( dconjg(geig1(1:ngp1,j1)) * qgg(1:ngp1,j,iwf))
       enddo
       enddo
       enddo

       deallocate(cphi1)
 1050 continue


c--- write result 
ctm
      do ispin = 1,nspin
         iti = iko_ixs(ispin)
         itf = iko_fxs(ispin)
         write(ifpsig(ispin)) iqbz
         write(ifpsig(ispin))
     &         ((psig(j1,j2,ispin),j1=iti,itf),j2=1,nwf)
      enddo

check write
      if (mod(iqbz,10).eq.1) then
         do ispin = 1,nspin
            iti = iko_ixs(ispin)
            itf = iko_fxs(ispin)
            do iwf = 1,nwf
               pgnorm = 0d0
               do j1 = iti,itf
                  pgnorm = pgnorm + abs(psig(j1,iwf,ispin))**2
c                  write(98,"(3i5,2f12.6)")iqbz,iwf,j1,
c     &                  abs(psig(j1,iwf,ispin))**2
               enddo
               write(*,*)ispin,iwf,pgnorm
            enddo
         enddo
      endif
      do ispin = 1,nspin
         iti = iko_ixs(ispin)
         itf = iko_fxs(ispin)
         do j1 = iti,itf
         do j2 = iti,itf
            sij = 0d0
            do iwf = 1,nwf
               sij = sij + dconjg(psig(j1,iwf,ispin))*psig(j2,iwf,ispin)
            enddo
            if (j1.eq.j2) then
               write(98,"(2i5,f12.6)")iqbz,j1,sij
            else
               write(97,"(3i5,f12.6)")iqbz,j1,j2,sij
            endif
         enddo
         enddo
      enddo

      deallocate(ngvecpf1, qgg, phig)

c end of qbz-loop
 1070 continue

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
     i             rr,nrofi,rmax,nphi,nc_max,
     d             nwf,nlmto,nrx,nl,nn,nbas,nspin,nphix,
     o             c1,c2)

      implicit real*8(a-h,o-z)
      real(8) :: phitotr(nrx,0:nl-1,nn,nbas,nspin),r0g(nphix,nwf),
     .           rr(nrx,nbas),rmax(nbas),
     .           c1(nphix,nwf,nspin),c2(nphix,nwf,nspin)
      integer(4) :: iphi(nphix,nwf),iphidot(nphix,nwf),
     .              n_indx(nlmto),l_indx(nlmto),m_indx(nlmto),
     .              ibas_indx(nlmto),nrofi(nbas),nphi(nwf),
     .              nc_max(0:nl-1,nbas)


      pi = 4d0*datan(1d0)

      c1 = 0d0
      c2 = 0d0

      do is  = 1,nspin
      do iwf = 1,nwf
      do j   = 1,nphi(iwf)

         ip = iphi(j,iwf)
         ibp = ibas_indx(ip)
         ilp = l_indx(ip)
         imp = m_indx(ip)
         inp = n_indx(ip) + nc_max(ilp,ibp)

         id = iphidot(j,iwf)
         ibd = ibas_indx(id)
         ild = l_indx(id)
         imd = m_indx(id)
         ind = n_indx(id) + nc_max(ild,ibd)

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
      subroutine phigmat(iphi,iphidot,c1,c2,phitoto,phitotr,
     i             rr,nrofi,rmax,aa,bb,r0g,ghead,
     i             n_indx,l_indx,m_indx,ibas_indx,nphi,nc_max,
     d             nwf,nlmto,nn,nl,nbas,nrx,nspin,nphix,
     o             phig)
c < phi | g>
c g = c1*phi + c2*phidot, which is conneted to the gaussian tail
c     at r=rmt

      implicit real*8 (a-h,o-z)

      real(8) :: phig(nlmto,nphix,nwf,nspin),
     &           phitoto(nrx,0:nl-1,nn,nbas,nspin),
     &           phitotr(nrx,0:nl-1,nn,nbas,nspin)
      real(8) :: rr(nrx,nbas),rmax(nbas),aa(nbas),bb(nbas),
     &           c1(nphix,nwf,nspin),c2(nphix,nwf,nspin),
     &           rphi(nrx),rhead(nrx),r0g(nphix,nwf)
      integer(4) :: iphi(nphix,nwf),iphidot(nphix,nwf),nrofi(nbas),
     &  n_indx(nlmto),l_indx(nlmto),m_indx(nlmto),ibas_indx(nlmto),
     &  nphi(nwf),nc_max(0:nl-1,nbas)
      logical :: ghead

      pi   = 4d0*datan(1d0)
      phig = 0d0

      do is  = 1,nspin
      do iwf = 1,nwf
      do j   = 1,nphi(iwf)
         
       if (ghead) then
         inlmp = iphi(j,iwf)
         ibasp = ibas_indx(inlmp)
         ilp   = l_indx(inlmp)
         imp   = m_indx(inlmp)
         inp   = n_indx(inlmp) + nc_max(ilp,ibasp)
c g(r) = A exp(-(r/r0)**2)
c A = (128 / (pi*r0))**(1/4)
         Ag = (128d0 / (pi*r0g(j,iwf)**6))**0.25d0
         rhead(1) = 0d0
         do ir = 2,nrofi(ibasp)
            rmt = rr(ir,ibasp)
            g  = Ag * dexp(-(rmt/r0g(j,iwf))**2)
            rhead(ir) = rmt * g
         enddo

       else ! ghead
c r * phi
         inlmp = iphi(j,iwf)
         ibasp = ibas_indx(inlmp)
         ilp   = l_indx(inlmp)
         imp   = m_indx(inlmp)
         inp   = n_indx(inlmp) + nc_max(ilp,ibasp)

c r * phidot
         inlmd = iphidot(j,iwf)
         ibasd = ibas_indx(inlmd)
         ild    = l_indx(inlmd)
         imd    = m_indx(inlmd)
         ind    = n_indx(inlmd) + nc_max(ild,ibasd)

         if (ibasp .ne. ibasd) stop 'phigmat: bas(phi) .ne. bas(phidot)'
         if (ilp .ne. ild) stop 'phigmat: l(phi) .ne. l(phidot)'
         if (imp .ne. imd) stop 'phigmat: m(phi) .ne. m(phidot)'

         rhead(1) = 0d0
         do ir = 2,nrofi(ibasp)
            rhead(ir) = c1(j,iwf,is) * phitotr(ir,ilp,inp,ibasp,is)
     &                + c2(j,iwf,is) * phitotr(ir,ild,ind,ibasd,is)
         enddo

       endif ! ghead

         do inlm = 1,nlmto
            ibas = ibas_indx(inlm)
            il   = l_indx(inlm)
            im   = m_indx(inlm)
            in   = n_indx(inlm) + nc_max(il,ibas)
            if (ibas .ne. ibasp) cycle
            if (il .ne. ilp) cycle
            if (im .ne. imp) cycle
            call gintxx(phitoto(1:nrofi(ibas),il,in,ibas,is),
     &                  rhead,aa(ibas),bb(ibas),nrofi(ibas),sum)
            phig(inlm,j,iwf,is) = sum
c            write(*,*)inlm,iwf,sum
         enddo

      enddo
c end of iwf-loop
      enddo
c end of is-loop
      enddo

c debug:
c      write(900,*)'   m    n    l   bas nlm   wf'
c      write(901,*)'   m    n    l   bas nlm   wf'
c      do iwf = 1,nwf
c      do inlm = 1,nlmto
c         ibas = ibas_indx(inlm)
c         in   = n_indx(inlm)
c         il   = l_indx(inlm)
c         im   = m_indx(inlm)
c         write(900,"(6i5,f12.6)")im,in,il,ibas,inlm,iwf,
c     &                            phig(inlm,1,iwf,1)
c         if (abs(phig(inlm,1,iwf,1)).gt.1.d-4)
c     &   write(901,"(6i5,f12.6)")im,in,il,ibas,inlm,iwf,
c     &                            phig(inlm,1,iwf,1)
c      enddo
c      enddo
c      write(901,*)'iwf, c1 c2'
c      do iwf=1,nwf
c         write(901,"(i5,2f12.6)")iwf,c1(1,iwf,1),c2(1,iwf,1)
c      enddo
c      stop 'debug: ok'

      return
      end
c-----------------------------------------------------------------------
