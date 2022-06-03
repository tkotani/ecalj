  !$$$      program h_uumatrix
  !$$$c-------------------------------------------------
  !$$$c  Calculate <u|u> matrix . u_kj(r) is the perodic part of eigencuntion.
  !$$$c This is a test routine. A bit detailed comment---so this code is a kind of manual
  !$$$c  to see the treatment of eigenfunctions.
  !$$$c
  !$$$c-------------------------------------------------
  !$$$      use m_readqg
  !$$$      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,readcphi,readgeig
  !$$$      use m_read_bzdata,ngrp2=>ngrp
  !$$$      use m_genallcf_v3
  !$$$      use keyvalue
  !$$$
  !$$$      implicit none
  !$$$      real(8):: q(3),  qgbin(3),qx(3)
  !$$$      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
  !$$$     & , ngvecpf1(:,:), ngvecpf2(:,:),
  !$$$     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
  !$$$      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:),
  !$$$     &phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:)
  !$$$      complex(8),parameter:: img=(0d0,1d0)
  !$$$c,nw,incwf,natom,nclass,ipos,igrp,
  !$$$c     & iinvg,nspin,nl,nn,nnv,nnc,
  !$$$c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions
  !$$$c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
  !$$$c     o                   izdummy,
  !$$$c     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
  !$$$c     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
  !$$$c     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
  !$$$c     o   iecore,ikonf,iicore,incore,nctot,             !core
  !$$$c     o   imagw_dummy,niw,idummy,
  !$$$      integer(4)
  !$$$     &   nw_input,
  !$$$     &   ifhbe,
  !$$$     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
  !$$$     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
  !$$$     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
  !$$$     &   ifepscond,nxx,ifvxcpout,ifgb0vec
  !$$$     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
  !$$$     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
  !$$$     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
  !$$$     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax  !,ifev1,ifev2
  !$$$     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
  !$$$     &   ,ngpmx
  !$$$
  !$$$      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini,
  !$$$     & dwry,dwh,omg_c,omg2,xxx
  !$$$      integer(4)::nwin, incwfin,  verbose
  !$$$      real(8)::efin
  !$$$      integer(4):: bzcase, mrecg,ifphi,
  !$$$     & nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc,
  !$$$     & ldim2,ixx,ngp1,ngp2,nq0it
  !$$$      real(8):: qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
  !$$$      real(8),parameter::  pi =     3.14159265358979323846d0
  !$$$      real(8),parameter::  fpi =    4d0*pi
  !$$$
  !$$$      logical:: qbzreg
  !$$$!-------------------------------------------------------------------------
  !$$$      integer(4),allocatable:: ncindx(:,:),
  !$$$     &           lcindx(:,:),
  !$$$     &           nrad(:),
  !$$$     &           nindx_r(:,:),
  !$$$     &           lindx_r(:,:),
  !$$$     &           nc_max(:,:),
  !$$$     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
  !$$$      real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:)
  !$$$     &                     ,phitotr(:,:,:,:,:),
  !$$$     &        bb(:),zz(:),rmax(:),cy(:),yl(:)
  !$$$
  !$$$
  !$$$      complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:)
  !$$$     & ,uum(:,:,:), ppovl(:,:)
  !$$$      complex(8):: ppj,phaseatom
  !$$$      real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
  !$$$      integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin
  !$$$     & ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx
  !$$$     &, iopen,ico,lxd,lx
  !$$$      real(8):: ylk
  !$$$#ifdef COMMONLL
  !$$$      integer(4):: ll(51**2)
  !$$$      common/llblock/ll
  !$$$#else
  !$$$      integer(4) ll
  !$$$#endif
  !$$$c-------------------------
  !$$$      call headver('h_uumatrix',0)
  !$$$c---  readin BZDATA. See gwsrc/rwbzdata.f
  !$$$c--------readin data set when you call read_BZDATA ---------------
  !$$$c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
  !$$$c     &   n_index_qbz
  !$$$c       integer(4):: n1,n2,n3
  !$$$c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_
  !$$$c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
  !$$$c     &    ,wibz(:),qbzw(:,:)
  !$$$c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
  !$$$c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
  !$$$c-----------------------------------------------------------------
  !$$$      call read_BZDATA()
  !$$$
  !$$$c--- Use regular mesh even for bzcase==2
  !$$$      if(bzcase()==2.and.qbzreg()) then
  !$$$        deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
  !$$$        do i=1,nqbz
  !$$$          qbz(:,i) = qbz(:,i) -deltaq/2d0
  !$$$          write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
  !$$$        enddo
  !$$$      endif
  !$$$      write(6,*)' ======== nqbz qbz  =',nqbz
  !$$$      write(6,*)' ======== nqibz ngrp=',nqibz,ngrp
  !$$$      write(6,*)  qbz
  !$$$      write(6,*)'============================'
  !$$$      print *
  !$$$
  !$$$C--- readin GWIN and LMTO, then allocate and set datas.
  !$$$      nwin = -999    !readin condition. Not readin NW file
  !$$$      incwfin= 0     !readin condition. use ForX0 for core in GWIN
  !$$$      efin =  -999d0 !readin condition. Not readin EFERMI
  !$$$      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
  !$$$Cstop2rx 2013.08.09 kino      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
  !$$$      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
  !$$$c---  These are allocated and setted by genallcf_v3
  !$$$c      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
  !$$$c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
  !$$$c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
  !$$$c      character(120):: symgrp
  !$$$c      character(6),allocatable :: clabl(:)
  !$$$c      integer(4),allocatable:: iclass(:)
  !$$$c     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
  !$$$c     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
  !$$$c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
  !$$$c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
  !$$$c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
  !$$$c     &    occv(:,:,:),unoccv(:,:,:)
  !$$$c     &   ,occc(:,:,:),unoccc(:,:,:),
  !$$$c     o    nocc(:,:,:),nunocc(:,:,:)
  !$$$c      real(8), allocatable::
  !$$$c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)
  !$$$
  !$$$!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
  !$$$Cstop2rx 2013.08.09 kino      if(nclass /= natom) stop ' nclass /= natom '
  !$$$      if(nclass /= natom) call rx( ' nclass /= natom ')
  !$$$
  !$$$c --- read dimensions of h,hb
  !$$$      ifhbe      = iopen('hbe.d',1,0,0)
  !$$$      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  !$$$Cstop2rx 2013.08.09 kino      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
  !$$$      if(nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
  !$$$Cstop2rx 2013.08.09 kino      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'
  !$$$      if(nqbz /=nqbzt ) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')
  !$$$
  !$$$c --- read by rdpp ; Radial integrals ppbrd and plane wave part
  !$$$      call getsrdpp2(nclass,nl,nxx)
  !$$$      call readngmx('QGpsi',ngpmx)
  !$$$      write(6,*)' ngpmx=',ngpmx
  !$$$
  !$$$c --- read radial functions PHIVC   (taken from hasfp0)
  !$$$      write(6,*)' Go to readining phivc'
  !$$$      ifphi  = iopen('PHIVC', 0,-1,0)     ! PHIV+PHIC augmentation wave and core
  !$$$      read(ifphi) nbas, nradmx, ncoremx,nrx
  !$$$Cstop2rx 2013.08.09 kino      if( nbas/=natom ) stop ' nbas(PHIVC) /= natom '
  !$$$      if( nbas/=natom ) call rx( ' nbas(PHIVC) /= natom ')
  !$$$      deallocate(ncore)
  !$$$      allocate(  ncindx(ncoremx,nbas),
  !$$$     &           lcindx(ncoremx,nbas),
  !$$$     &           nrad(nbas),
  !$$$     &           nindx_r(1:nradmx,1:nbas),
  !$$$     &           lindx_r(1:nradmx,1:nbas),
  !$$$     &        aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) ,
  !$$$     &        phitoto(nrx,0:nl-1,nn,nbas,nspin),
  !$$$     &        phitotr(nrx,0:nl-1,nn,nbas,nspin),
  !$$$     &        nc_max(0:nl-1,nbas),ncore(nbas),rmax(nbas) )
  !$$$      write(6,*)' end of allocation'
  !$$$      read(ifphi) nrad(1:nbas)
  !$$$      read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
  !$$$      nc_max=0
  !$$$      do ibas=1,nbas
  !$$$        ic = ibas
  !$$$        write(6,*)' --- read PHIVC of ibas nrad=',ibas,nrad(ic)
  !$$$        read(ifphi) ncore(ic), ncoremx                            !core
  !$$$        read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
  !$$$        write(6,*)' xxx0'
  !$$$        read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
  !$$$
  !$$$        write(6,*) 'ic icx=',ic,icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
  !$$$        if(ic/=icx) then
  !$$$Cstop2rx 2013.08.09 kino          stop ' h_uu: ic/=icx'
  !$$$          call rx( ' h_uu: ic/=icx')
  !$$$        endif
  !$$$        write(6,*)' xxx1 ncoremx ncore(ic)=',ncoremx,ncore(ic)
  !$$$        read(ifphi) rr(1:nrofi(ic),ic)
  !$$$        write(6,*)' xxx2 ncoremx ncore(ic)=',ncoremx,ncore(ic)
  !$$$
  !$$$        write(6,*)' xxx2 nspin=',nspin
  !$$$        rmax(ic) = rr(nrofi(ic),ic)
  !$$$        do isp = 1, nspin
  !$$$          write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
  !$$$          do ico = 1, ncore(ic) !core
  !$$$            l =  lcindx(ico,ic)
  !$$$            n =  ncindx(ico,ic)
  !$$$            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
  !$$$            phitotr(1:nrofi(ic),l,n, ic,isp)=              !core raw= core orthgonal
  !$$$     &      phitoto(1:nrofi(ic),l,n, ic,isp)               !
  !$$$            if(n>nc_max(l,ic)) nc_max(l,ic)=n
  !$$$            write(6,*)' sss1c=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
  !$$$          enddo
  !$$$          do irad = 1, nrad(ic)   !valence
  !$$$            l = lindx_r (irad,ic)
  !$$$            n = nindx_r (irad,ic) + nc_max(l,ic)
  !$$$            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
  !$$$            read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
  !$$$            write(6,*)' sss1=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
  !$$$            write(6,*)' sss2=',sum(abs(phitotr(1:nrofi(ic),l,n, ic,isp)))
  !$$$          enddo
  !$$$        enddo
  !$$$      enddo
  !$$$
  !$$$c--- cg coefficient.  y = cg y y ; y is the real spherical harmonics
  !$$$      ngrpx=1
  !$$$      allocate( cg(nl**2,nl**2,(2*nl-1)**2), symope(3,3) )
  !$$$      symope(1:3,1) = (/1d0,0d0,0d0/)
  !$$$      symope(1:3,2) = (/0d0,1d0,0d0/)
  !$$$      symope(1:3,3) = (/0d0,0d0,1d0/)
  !$$$      cg = 0d0 !for sanity check
  !$$$      call rotcg(nl-1,symope,ngrpx,cg)
  !$$$
  !$$$c --- initiallization to get eigenfunctions
  !$$$      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
  !$$$      call init_readeigen2(mrecb,nlmto,mrecg)
  !$$$      call readngmx('QGpsi',ngpmx)
  !$$$      allocate( geig1(ngpmx,nband),geig2(ngpmx,nband))
  !$$$      write(6,*) 'end of initialization'
  !$$$
  !$$$c --- Readin nlam index
  !$$$      ifoc = iopen('@MNLA_CPHI',1,0,0)
  !$$$      ldim2 = nlmto
  !$$$      read(ifoc,*)
  !$$$      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
  !$$$      do ix =1,ldim2
  !$$$        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
  !$$$Cstop2rx 2013.08.09 kino        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
  !$$$        if(ixx/=ix) call rx( 'failed to readin @MNLA_CPHI')
  !$$$      enddo
  !$$$
  !$$$c ---  q near zero
  !$$$      write(6,*) 'reading QOP'
  !$$$      open (101,file='Q0P')
  !$$$      read (101,"(i5)") nq0i
  !$$$!      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
  !$$$      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
  !$$$      nq0it = nq0i
  !$$$      allocate( q0i(1:3,1:nq0i) ) !wqt(1:nq0i),
  !$$$!      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
  !$$$      nq0ix = nq0i
  !$$$      do i=1,nq0i
  !$$$        read (101,* ) xxx,q0i(1:3,i)
  !$$$        if(xxx==0d0 ) nq0ix = i-1
  !$$$      enddo
  !$$$      nq0i = nq0ix ! New nq0i July 2001
  !$$$      write(6,*) ' Used k number in Q0P =', nq0i
  !$$$      write(6,"(i3, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
  !$$$      close(101)
  !$$$
  !$$$
  !$$$c======================================================================
  !$$$c --- Set q1(j1range) q2(j2range)
  !$$$c======================================================================
  !$$$! Note that the true q when we generate eigenfunctions are q1x and q2x.
  !$$$! q1-q1x should be a G vector.
  !$$$! So you may need to take into account the phase shift to <u|u> vectors.
  !$$$!
  !$$$! --- I inserted checkagree to make sure that q1=q1x and q2=q2x ...
  !$$$      q1 = qbz(:, 1) ;  j1min=1;   j1max=8
  !$$$      q2 = qbz(:, 2) ;  j2min=1;   j2max=8
  !$$$cc      q1 = qbz(:, 12)           ;  j1min=1;   j1max=8
  !$$$cc      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
  !$$$c      q1 = qbz(:, 5);  j1min=1;   j1max=8
  !$$$c      q2 = qbz(:, 5);  j2min=1;   j2max=8
  !$$$c      q1 = qbz(:, 12)+q0i(1:3,1);  j1min=1;   j1max=8
  !$$$c      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
  !$$$c======================================================================
  !$$$      allocate( uum(j1min:j1max,j2min:j2max,nspin) )
  !$$$
  !$$$C --- ppovl= <P_{q1+G1}|P_{q2+G2}>
  !$$$      call readqg0('QGpsi',q1,ginv,  q1x, ngp1)
  !$$$      call readqg0('QGpsi',q2,ginv,  q2x, ngp2)
  !$$$      call checkagree(q1,q1x,' q1 ne q1x') ! make sure q1=q1x
  !$$$      call checkagree(q2,q2x,' q2 ne q2x')
  !$$$      write(6,"('q1 =',3f9.4,3x,3f9.4)") q1
  !$$$      write(6,"('q2 =',3f9.4,3x,3f9.4)") q2
  !$$$      allocate( ngvecpf1(3,ngp1), ngvecpf2(3,ngp2), ppovl(ngp1,ngp2) )
  !$$$      call readqg('QGpsi',q1,ginv, q1x, ngp1, ngvecpf1)
  !$$$      call readqg('QGpsi',q2,ginv, q2x, ngp2, ngvecpf2)
  !$$$      call checkagree(q1,q1x,' q1 ne q1x xxx2') !make sure q1 =q1x
  !$$$      call checkagree(q2,q2x,' q2 ne q2x xxx2')
  !$$$      write(6,*)' ngp1,ngp2 sum check=',ngp1,ngp2,sum(abs(ngvecpf1)),sum(abs(ngvecpf2))
  !$$$      call mkppovl2(alat,plat,qbas, !22April2004
  !$$$     &    ngp1, ngvecpf1,
  !$$$     &    ngp2, ngvecpf2,
  !$$$     &    nbas, rmax, pos,
  !$$$     o    ppovl)
  !$$$      write(6,*)' end of mkppovl2'
  !$$$
  !$$$c ... lxx and allocations
  !$$$      lxx=2*(nl-1)
  !$$$      allocate( ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nspin,nbas),
  !$$$     &   rprodx(nrx,0:lxx),
  !$$$     &   phij(0:lxx),psij(0:lxx),rphiphi(nrx))
  !$$$
  !$$$c ... dq
  !$$$      dq = q1x-q2x
  !$$$      if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
  !$$$
  !$$$      absdq = sqrt(sum(dq**2))
  !$$$      absqg2 = (2*pi/alat)**2 *sum(dq**2)
  !$$$      absqg =sqrt(absqg2)
  !$$$
  !$$$c ... YL(dq)
  !$$$      allocate(cy((lxx+1)**2),yl((lxx+1)**2))
  !$$$      call sylmnc(cy,lxx)
  !$$$      call sylm(dq/absdq,yl,lxx,r2s) !spherical factor Y(dq)
  !$$$
  !$$$C --- radial integral  ppbrd = <phi phi j_l>
  !$$$      ppbrd=0d0
  !$$$      do 900 ibas = 1,nbas
  !$$$        ic = ibas
  !$$$        write(6,"(' nindx=',10i3)") nindx(1:nl,ic)
  !$$$        write(6,*)' radial integral ibas=',ibas
  !$$$        do ir =1,nrofi(ic)
  !$$$          call bessl(absqg2*rr(ir,ibas)**2,lxx,phij,psij)
  !$$$c  phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
  !$$$          do lx = 0, lxx
  !$$$            rprodx(ir,lx) = rr(ir,ibas)* phij(lx)* (absqg*rr(ir,ibas))**lx
  !$$$            ! = r \times j_l(|dq|r)  !bessel function
  !$$$          enddo
  !$$$ccccccccccccccccccccccc
  !$$$c          write(1100,"(10d13.6)")rr(ir,ibas),rprodx(ir,0:lxx) !,phij(0:lxx)
  !$$$ccccccccccccccccccccccc
  !$$$        enddo
  !$$$        do 125 isp = 1,nspin
  !$$$          do 25 l1 = 0, nl-1
  !$$$          do 25 n1 = 1, nindx(l1+1,ic)
  !$$$          do 25 l2 = 0, nl-1
  !$$$          do 25 n2 = 1, nindx(l2+1,ic)
  !$$$            rphiphi(1)       = 0d0
  !$$$            rphiphi(2:nrofi(ic)) = phitoto(2:nrofi(ic),l1,n1,ic,isp)
  !$$$     &                          *phitoto(2:nrofi(ic),l2,n2,ic,isp)/rr(2:,ic) ! phi = u = r \phi
  !$$$          do 25 lx = 0, 2*(nl-1)
  !$$$            if(lx <abs(l1-l2) .or. l1+l2<lx) cycle
  !$$$            call gintxx( rprodx(1,lx), rphiphi,aa(ic),bb(ic),nrofi(ic),
  !$$$     &        ppbrd(l1, n1,l2, n2, lx, isp,ibas) )
  !$$$c          if(l1==l2.and.n1==n2.and.lx==0)
  !$$$c         write(6,*) ' ppbrd=',l1,n1,ppbrd(l1, n1,l2, n2, lx, isp,ibas)
  !$$$   25     continue
  !$$$ 125    continue
  !$$$ 900    continue
  !$$$
  !$$$C --- Calcuate <u{q1x j1} | u_{q2x j2}>
  !$$$c              = < exp(i(q1x-q2x)r) psi^*{q1x j1} psi_{q2x j2} >
  !$$$c ... MT part
  !$$$cr   ldim2 = nlmto
  !$$$cr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
  !$$$cr   l_indx   (1:ldim2) : l index
  !$$$cr   ibas_indx(1:ldim2) : ibas index.
  !$$$        uum = 0d0
  !$$$        do 1050 ispin=1,nspin
  !$$$          allocate(cphi1 (nlmto,nband),cphi2(nlmto,nband) )
  !$$$          call readcphi(q1, nlmto, ispin, quu, cphi1) !quu is true q, q1-q can be G vectors. But need to check it again!
  !$$$                                                      !Need to examine readcphi again.
  !$$$          call checkagree(q1,q1x,' q1 ne quu') !check for safe quu is true q
  !$$$          call readcphi(q2, nlmto, ispin, quu, cphi2)
  !$$$          call checkagree(q2,quu,' q2 ne quu')
  !$$$
  !$$$          do 1020 ia1 = 1,nlmto
  !$$$            ibas1= ibas_indx(ia1)
  !$$$            l1   = l_indx    (ia1)
  !$$$            m1   = m_indx    (ia1)
  !$$$            n1   = n_indx    (ia1) + nc_max(l1,ibas1)
  !$$$            lm1  = l1**2+l1+1  + m1
  !$$$            do 1010 ia2 = 1,nlmto
  !$$$              ibas2 = ibas_indx(ia2)
  !$$$              if(ibas2/=ibas1) cycle
  !$$$              phaseatom = exp( img* 2d0*pi*sum(dq*pos(:,ibas1)) )
  !$$$              m2   = m_indx    (ia2)
  !$$$              l2   = l_indx    (ia2)
  !$$$              n2   = n_indx    (ia2) + nc_max(l2,ibas2)
  !$$$              lm2= l2**2 +l2+1 + m2
  !$$$cccccccccccccccccccccccccccc
  !$$$c Norm check test.
  !$$$c          do j1= j1min,j1max
  !$$$c          do j2= j2min,j2max
  !$$$c            if(ia1==ia2) uum(j1,j2,ispin) = uum(j1,j2,ispin)
  !$$$c     &        + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2)
  !$$$c          enddo
  !$$$c          enddo
  !$$$ccccccccccccccccccccccccccccc
  !$$$              do lm3= (l1-l2)**2+1, (l1+l2+1)**2 ! l3 can take |l1-l2|,...l1+l2
  !$$$                l3 = ll(lm3)
  !$$$ccccccccccccccccccccccccccccccc
  !$$$c          ylk=0d0;   if(l3==0) ylk=1d0/sqrt(4*pi) !Y_00 only test
  !$$$ccccccccccccccccccccccccccccccc
  !$$$                ylk= cy(lm3)*yl(lm3)
  !$$$                ppj = ppbrd(l1,n1,l2,n2,l3,ispin,ibas1) *cg(lm1,lm2, lm3)
  !$$$     &          * fpi* img**l3* phaseatom * ylk
  !$$$
  !$$$! cg(lm1,lm2,lm3)= \int Y_lm3(\hat(r)) Y_lm2(\hat(r)) Y_lm1(\hat(r)) \frac{d \Omega}{4\pi}
  !$$$! This is based on inverse expansion. See Rose.Eq.3.8.
  !$$$                do j1= j1min,j1max
  !$$$                  do j2= j2min,j2max
  !$$$                    uum(j1,j2,ispin) =
  !$$$     &      uum(j1,j2,ispin) + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2) * ppj
  !$$$                  enddo
  !$$$                enddo
  !$$$              enddo
  !$$$ 1010       continue
  !$$$ 1020     continue
  !$$$
  !$$$c ... Interstitial Plane Wave part
  !$$$          call readgeig(q1, ngpmx, ispin, quu, geig1)
  !$$$          call checkagree(q1,quu,' q1 ne quu eig')
  !$$$          call readgeig(q2, ngpmx, ispin, quu, geig2)
  !$$$          call checkagree(q2,quu,' q1 ne quu eig')
  !$$$          do j1=j1min,j1max
  !$$$            do j2=j2min,j2max
  !$$$              uum(j1,j2,ispin)= uum(j1,j2,ispin) +
  !$$$     &    sum( dconjg(geig1(1:ngp1,j1)) * matmul(ppovl,geig2(1:ngp2,j2)) )
  !$$$            enddo
  !$$$          enddo
  !$$$        deallocate(cphi1,cphi2)
  !$$$ 1050   continue
  !$$$
  !$$$c--- write resutlt
  !$$$        write(6,*)' ============ result --- diagonal --- =============='
  !$$$        do ispin = 1,nspin
  !$$$          do j1=j1min,j1max
  !$$$            do j2=j2min,j2max
  !$$$              if(j1==j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)")
  !$$$     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
  !$$$            enddo
  !$$$          enddo
  !$$$          write(6,*)'--- off diagonal ----------------------'
  !$$$          do j1=j1min,j1max
  !$$$            do j2=j2min,j2max
  !$$$              if(j1/=j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)")
  !$$$     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
  !$$$            enddo
  !$$$          enddo
  !$$$        enddo
  !$$$        write(6,*) ' ====== end ========================================'
  !$$$c      stop ' ====== end ========================================'
  !$$$        end
  !$$$
  !$$$        subroutine checkagree(a,b,char)
  !$$$        real(8):: a(3),b(3)
  !$$$        character*(*) :: char
  !$$$        if(sum(abs(a-b))>1d-6) then
  !$$$          write(6,*)' Error in checkagree:',char
  !$$$Cstop2rx 2013.08.09 kino          stop ' Error in checkagree:'
  !$$$          call rx( ' Error in checkagree:')
  !$$$        endif
  !$$$        end
  