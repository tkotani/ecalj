c$$$      program h_uumatrix
c$$$c-------------------------------------------------
c$$$c  Calculate <u|u> matrix . u_kj(r) is the perodic part of eigencuntion.
c$$$c This is a test routine. A bit detailed comment---so this code is a kind of manual
c$$$c  to see the treatment of eigenfunctions.
c$$$c
c$$$c-------------------------------------------------
c$$$      use m_readqg
c$$$      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,readcphi,readgeig
c$$$      use m_read_bzdata,ngrp2=>ngrp
c$$$      use m_genallcf_v3
c$$$      use keyvalue
c$$$
c$$$      implicit none
c$$$      real(8):: q(3),  qgbin(3),qx(3)
c$$$      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
c$$$     & , ngvecpf1(:,:), ngvecpf2(:,:), 
c$$$     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
c$$$      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:),
c$$$     &phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:)
c$$$      complex(8),parameter:: img=(0d0,1d0)
c$$$c,nw,incwf,natom,nclass,ipos,igrp,
c$$$c     & iinvg,nspin,nl,nn,nnv,nnc,
c$$$c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions
c$$$c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
c$$$c     o                   izdummy,
c$$$c     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
c$$$c     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
c$$$c     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
c$$$c     o   iecore,ikonf,iicore,incore,nctot,             !core
c$$$c     o   imagw_dummy,niw,idummy,
c$$$      integer(4)
c$$$     &   nw_input,
c$$$     &   ifhbe,
c$$$     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
c$$$     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
c$$$     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
c$$$     &   ifepscond,nxx,ifvxcpout,ifgb0vec
c$$$     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
c$$$     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
c$$$     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
c$$$     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax  !,ifev1,ifev2
c$$$     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
c$$$     &   ,ngpmx
c$$$
c$$$      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini,
c$$$     & dwry,dwh,omg_c,omg2,xxx
c$$$      integer(4)::nwin, incwfin,  verbose
c$$$      real(8)::efin
c$$$      integer(4):: bzcase, mrecg,ifphi,
c$$$     & nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc,
c$$$     & ldim2,ixx,ngp1,ngp2,nq0it
c$$$      real(8):: qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
c$$$      real(8),parameter::  pi =     3.14159265358979323846d0
c$$$      real(8),parameter::  fpi =    4d0*pi
c$$$
c$$$      logical:: qbzreg
c$$$!-------------------------------------------------------------------------
c$$$      integer(4),allocatable:: ncindx(:,:),
c$$$     &           lcindx(:,:),
c$$$     &           nrad(:),
c$$$     &           nindx_r(:,:),
c$$$     &           lindx_r(:,:),
c$$$     &           nc_max(:,:),
c$$$     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
c$$$      real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:)
c$$$     &                     ,phitotr(:,:,:,:,:),
c$$$     &        bb(:),zz(:),rmax(:),cy(:),yl(:)
c$$$
c$$$
c$$$      complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:)
c$$$     & ,uum(:,:,:), ppovl(:,:)
c$$$      complex(8):: ppj,phaseatom
c$$$      real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
c$$$      integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin
c$$$     & ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx
c$$$     &, iopen,ico,lxd,lx
c$$$      real(8):: ylk
c$$$#ifdef COMMONLL
c$$$      integer(4):: ll(51**2)
c$$$      common/llblock/ll
c$$$#else
c$$$      integer(4) ll
c$$$#endif
c$$$c-------------------------
c$$$      call headver('h_uumatrix',0)
c$$$c---  readin BZDATA. See gwsrc/rwbzdata.f
c$$$c--------readin data set when you call read_BZDATA ---------------
c$$$c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c$$$c     &   n_index_qbz
c$$$c       integer(4):: n1,n2,n3
c$$$c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_
c$$$c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c$$$c     &    ,wibz(:),qbzw(:,:)
c$$$c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c$$$c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
c$$$c-----------------------------------------------------------------
c$$$      call read_BZDATA()
c$$$
c$$$c--- Use regular mesh even for bzcase==2
c$$$      if(bzcase()==2.and.qbzreg()) then
c$$$        deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
c$$$        do i=1,nqbz
c$$$          qbz(:,i) = qbz(:,i) -deltaq/2d0
c$$$          write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
c$$$        enddo
c$$$      endif
c$$$      write(6,*)' ======== nqbz qbz  =',nqbz
c$$$      write(6,*)' ======== nqibz ngrp=',nqibz,ngrp
c$$$      write(6,*)  qbz
c$$$      write(6,*)'============================'
c$$$      print *
c$$$
c$$$C--- readin GWIN and LMTO, then allocate and set datas.
c$$$      nwin = -999    !readin condition. Not readin NW file
c$$$      incwfin= 0     !readin condition. use ForX0 for core in GWIN
c$$$      efin =  -999d0 !readin condition. Not readin EFERMI
c$$$      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
c$$$Cstop2rx 2013.08.09 kino      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c$$$      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
c$$$c---  These are allocated and setted by genallcf_v3
c$$$c      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c$$$c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
c$$$c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
c$$$c      character(120):: symgrp
c$$$c      character(6),allocatable :: clabl(:)
c$$$c      integer(4),allocatable:: iclass(:)
c$$$c     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
c$$$c     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c$$$c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c$$$c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c$$$c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c$$$c     &    occv(:,:,:),unoccv(:,:,:)
c$$$c     &   ,occc(:,:,:),unoccc(:,:,:),
c$$$c     o    nocc(:,:,:),nunocc(:,:,:)
c$$$c      real(8), allocatable::
c$$$c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)
c$$$
c$$$!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
c$$$Cstop2rx 2013.08.09 kino      if(nclass /= natom) stop ' nclass /= natom '
c$$$      if(nclass /= natom) call rx( ' nclass /= natom ')
c$$$
c$$$c --- read dimensions of h,hb
c$$$      ifhbe      = iopen('hbe.d',1,0,0)
c$$$      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
c$$$Cstop2rx 2013.08.09 kino      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
c$$$      if(nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
c$$$Cstop2rx 2013.08.09 kino      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'
c$$$      if(nqbz /=nqbzt ) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')
c$$$
c$$$c --- read by rdpp ; Radial integrals ppbrd and plane wave part
c$$$      call getsrdpp2(nclass,nl,nxx)
c$$$      call readngmx('QGpsi',ngpmx)
c$$$      write(6,*)' ngpmx=',ngpmx
c$$$
c$$$c --- read radial functions PHIVC   (taken from hasfp0)
c$$$      write(6,*)' Go to readining phivc'
c$$$      ifphi  = iopen('PHIVC', 0,-1,0)     ! PHIV+PHIC augmentation wave and core
c$$$      read(ifphi) nbas, nradmx, ncoremx,nrx
c$$$Cstop2rx 2013.08.09 kino      if( nbas/=natom ) stop ' nbas(PHIVC) /= natom '
c$$$      if( nbas/=natom ) call rx( ' nbas(PHIVC) /= natom ')
c$$$      deallocate(ncore)
c$$$      allocate(  ncindx(ncoremx,nbas),
c$$$     &           lcindx(ncoremx,nbas),
c$$$     &           nrad(nbas),
c$$$     &           nindx_r(1:nradmx,1:nbas),
c$$$     &           lindx_r(1:nradmx,1:nbas),
c$$$     &        aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) ,
c$$$     &        phitoto(nrx,0:nl-1,nn,nbas,nspin),
c$$$     &        phitotr(nrx,0:nl-1,nn,nbas,nspin),
c$$$     &        nc_max(0:nl-1,nbas),ncore(nbas),rmax(nbas) )
c$$$      write(6,*)' end of allocation'
c$$$      read(ifphi) nrad(1:nbas)
c$$$      read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
c$$$      nc_max=0
c$$$      do ibas=1,nbas
c$$$        ic = ibas
c$$$        write(6,*)' --- read PHIVC of ibas nrad=',ibas,nrad(ic)
c$$$        read(ifphi) ncore(ic), ncoremx                            !core
c$$$        read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
c$$$        write(6,*)' xxx0'
c$$$        read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
c$$$
c$$$        write(6,*) 'ic icx=',ic,icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
c$$$        if(ic/=icx) then
c$$$Cstop2rx 2013.08.09 kino          stop ' h_uu: ic/=icx'
c$$$          call rx( ' h_uu: ic/=icx')
c$$$        endif
c$$$        write(6,*)' xxx1 ncoremx ncore(ic)=',ncoremx,ncore(ic)
c$$$        read(ifphi) rr(1:nrofi(ic),ic)
c$$$        write(6,*)' xxx2 ncoremx ncore(ic)=',ncoremx,ncore(ic)
c$$$
c$$$        write(6,*)' xxx2 nspin=',nspin
c$$$        rmax(ic) = rr(nrofi(ic),ic)
c$$$        do isp = 1, nspin
c$$$          write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
c$$$          do ico = 1, ncore(ic) !core
c$$$            l =  lcindx(ico,ic)
c$$$            n =  ncindx(ico,ic)
c$$$            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
c$$$            phitotr(1:nrofi(ic),l,n, ic,isp)=              !core raw= core orthgonal
c$$$     &      phitoto(1:nrofi(ic),l,n, ic,isp)               !
c$$$            if(n>nc_max(l,ic)) nc_max(l,ic)=n
c$$$            write(6,*)' sss1c=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
c$$$          enddo
c$$$          do irad = 1, nrad(ic)   !valence
c$$$            l = lindx_r (irad,ic)
c$$$            n = nindx_r (irad,ic) + nc_max(l,ic)
c$$$            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
c$$$            read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
c$$$            write(6,*)' sss1=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
c$$$            write(6,*)' sss2=',sum(abs(phitotr(1:nrofi(ic),l,n, ic,isp)))
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c--- cg coefficient.  y = cg y y ; y is the real spherical harmonics
c$$$      ngrpx=1
c$$$      allocate( cg(nl**2,nl**2,(2*nl-1)**2), symope(3,3) )
c$$$      symope(1:3,1) = (/1d0,0d0,0d0/)
c$$$      symope(1:3,2) = (/0d0,1d0,0d0/)
c$$$      symope(1:3,3) = (/0d0,0d0,1d0/)
c$$$      cg = 0d0 !for sanity check
c$$$      call rotcg(nl-1,symope,ngrpx,cg)
c$$$
c$$$c --- initiallization to get eigenfunctions
c$$$      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
c$$$      call init_readeigen2(mrecb,nlmto,mrecg)
c$$$      call readngmx('QGpsi',ngpmx)
c$$$      allocate( geig1(ngpmx,nband),geig2(ngpmx,nband))
c$$$      write(6,*) 'end of initialization'
c$$$
c$$$c --- Readin nlam index
c$$$      ifoc = iopen('@MNLA_CPHI',1,0,0)
c$$$      ldim2 = nlmto
c$$$      read(ifoc,*)
c$$$      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
c$$$      do ix =1,ldim2
c$$$        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
c$$$Cstop2rx 2013.08.09 kino        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
c$$$        if(ixx/=ix) call rx( 'failed to readin @MNLA_CPHI')
c$$$      enddo
c$$$
c$$$c ---  q near zero
c$$$      write(6,*) 'reading QOP'
c$$$      open (101,file='Q0P')
c$$$      read (101,"(i5)") nq0i
c$$$!      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
c$$$      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
c$$$      nq0it = nq0i
c$$$      allocate( q0i(1:3,1:nq0i) ) !wqt(1:nq0i),
c$$$!      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
c$$$      nq0ix = nq0i
c$$$      do i=1,nq0i
c$$$        read (101,* ) xxx,q0i(1:3,i)
c$$$        if(xxx==0d0 ) nq0ix = i-1
c$$$      enddo
c$$$      nq0i = nq0ix ! New nq0i July 2001
c$$$      write(6,*) ' Used k number in Q0P =', nq0i
c$$$      write(6,"(i3, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
c$$$      close(101)
c$$$
c$$$
c$$$c======================================================================
c$$$c --- Set q1(j1range) q2(j2range)
c$$$c======================================================================
c$$$! Note that the true q when we generate eigenfunctions are q1x and q2x.
c$$$! q1-q1x should be a G vector.
c$$$! So you may need to take into account the phase shift to <u|u> vectors.
c$$$!
c$$$! --- I inserted checkagree to make sure that q1=q1x and q2=q2x ...
c$$$      q1 = qbz(:, 1) ;  j1min=1;   j1max=8
c$$$      q2 = qbz(:, 2) ;  j2min=1;   j2max=8
c$$$cc      q1 = qbz(:, 12)           ;  j1min=1;   j1max=8
c$$$cc      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
c$$$c      q1 = qbz(:, 5);  j1min=1;   j1max=8
c$$$c      q2 = qbz(:, 5);  j2min=1;   j2max=8
c$$$c      q1 = qbz(:, 12)+q0i(1:3,1);  j1min=1;   j1max=8
c$$$c      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
c$$$c======================================================================
c$$$      allocate( uum(j1min:j1max,j2min:j2max,nspin) )
c$$$
c$$$C --- ppovl= <P_{q1+G1}|P_{q2+G2}>
c$$$      call readqg0('QGpsi',q1,ginv,  q1x, ngp1)
c$$$      call readqg0('QGpsi',q2,ginv,  q2x, ngp2)
c$$$      call checkagree(q1,q1x,' q1 ne q1x') ! make sure q1=q1x
c$$$      call checkagree(q2,q2x,' q2 ne q2x') 
c$$$      write(6,"('q1 =',3f9.4,3x,3f9.4)") q1
c$$$      write(6,"('q2 =',3f9.4,3x,3f9.4)") q2
c$$$      allocate( ngvecpf1(3,ngp1), ngvecpf2(3,ngp2), ppovl(ngp1,ngp2) )
c$$$      call readqg('QGpsi',q1,ginv, q1x, ngp1, ngvecpf1)
c$$$      call readqg('QGpsi',q2,ginv, q2x, ngp2, ngvecpf2)
c$$$      call checkagree(q1,q1x,' q1 ne q1x xxx2') !make sure q1 =q1x
c$$$      call checkagree(q2,q2x,' q2 ne q2x xxx2')
c$$$      write(6,*)' ngp1,ngp2 sum check=',ngp1,ngp2,sum(abs(ngvecpf1)),sum(abs(ngvecpf2))
c$$$      call mkppovl2(alat,plat,qbas, !22April2004 
c$$$     &    ngp1, ngvecpf1, 
c$$$     &    ngp2, ngvecpf2, 
c$$$     &    nbas, rmax, pos,  
c$$$     o    ppovl)
c$$$      write(6,*)' end of mkppovl2'
c$$$
c$$$c ... lxx and allocations
c$$$      lxx=2*(nl-1)
c$$$      allocate( ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nspin,nbas),
c$$$     &   rprodx(nrx,0:lxx),
c$$$     &   phij(0:lxx),psij(0:lxx),rphiphi(nrx))
c$$$
c$$$c ... dq
c$$$      dq = q1x-q2x
c$$$      if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
c$$$
c$$$      absdq = sqrt(sum(dq**2))
c$$$      absqg2 = (2*pi/alat)**2 *sum(dq**2)
c$$$      absqg =sqrt(absqg2)
c$$$
c$$$c ... YL(dq)
c$$$      allocate(cy((lxx+1)**2),yl((lxx+1)**2))
c$$$      call sylmnc(cy,lxx)
c$$$      call sylm(dq/absdq,yl,lxx,r2s) !spherical factor Y(dq)
c$$$
c$$$C --- radial integral  ppbrd = <phi phi j_l>
c$$$      ppbrd=0d0
c$$$      do 900 ibas = 1,nbas
c$$$        ic = ibas
c$$$        write(6,"(' nindx=',10i3)") nindx(1:nl,ic)
c$$$        write(6,*)' radial integral ibas=',ibas
c$$$        do ir =1,nrofi(ic)
c$$$          call bessl(absqg2*rr(ir,ibas)**2,lxx,phij,psij)
c$$$c  phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
c$$$          do lx = 0, lxx
c$$$            rprodx(ir,lx) = rr(ir,ibas)* phij(lx)* (absqg*rr(ir,ibas))**lx
c$$$            ! = r \times j_l(|dq|r)  !bessel function
c$$$          enddo
c$$$ccccccccccccccccccccccc
c$$$c          write(1100,"(10d13.6)")rr(ir,ibas),rprodx(ir,0:lxx) !,phij(0:lxx)
c$$$ccccccccccccccccccccccc
c$$$        enddo
c$$$        do 125 isp = 1,nspin
c$$$          do 25 l1 = 0, nl-1
c$$$          do 25 n1 = 1, nindx(l1+1,ic)
c$$$          do 25 l2 = 0, nl-1
c$$$          do 25 n2 = 1, nindx(l2+1,ic)
c$$$            rphiphi(1)       = 0d0
c$$$            rphiphi(2:nrofi(ic)) = phitoto(2:nrofi(ic),l1,n1,ic,isp)
c$$$     &                          *phitoto(2:nrofi(ic),l2,n2,ic,isp)/rr(2:,ic) ! phi = u = r \phi
c$$$          do 25 lx = 0, 2*(nl-1)
c$$$            if(lx <abs(l1-l2) .or. l1+l2<lx) cycle
c$$$            call gintxx( rprodx(1,lx), rphiphi,aa(ic),bb(ic),nrofi(ic),
c$$$     &        ppbrd(l1, n1,l2, n2, lx, isp,ibas) )
c$$$c          if(l1==l2.and.n1==n2.and.lx==0)
c$$$c         write(6,*) ' ppbrd=',l1,n1,ppbrd(l1, n1,l2, n2, lx, isp,ibas)
c$$$   25     continue
c$$$ 125    continue
c$$$ 900    continue
c$$$
c$$$C --- Calcuate <u{q1x j1} | u_{q2x j2}>
c$$$c              = < exp(i(q1x-q2x)r) psi^*{q1x j1} psi_{q2x j2} >
c$$$c ... MT part
c$$$cr   ldim2 = nlmto
c$$$cr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
c$$$cr   l_indx   (1:ldim2) : l index
c$$$cr   ibas_indx(1:ldim2) : ibas index.
c$$$        uum = 0d0
c$$$        do 1050 ispin=1,nspin
c$$$          allocate(cphi1 (nlmto,nband),cphi2(nlmto,nband) )
c$$$          call readcphi(q1, nlmto, ispin, quu, cphi1) !quu is true q, q1-q can be G vectors. But need to check it again!
c$$$                                                      !Need to examine readcphi again.  
c$$$          call checkagree(q1,q1x,' q1 ne quu') !check for safe quu is true q 
c$$$          call readcphi(q2, nlmto, ispin, quu, cphi2)
c$$$          call checkagree(q2,quu,' q2 ne quu')
c$$$
c$$$          do 1020 ia1 = 1,nlmto
c$$$            ibas1= ibas_indx(ia1)
c$$$            l1   = l_indx    (ia1)
c$$$            m1   = m_indx    (ia1)
c$$$            n1   = n_indx    (ia1) + nc_max(l1,ibas1)
c$$$            lm1  = l1**2+l1+1  + m1
c$$$            do 1010 ia2 = 1,nlmto
c$$$              ibas2 = ibas_indx(ia2)
c$$$              if(ibas2/=ibas1) cycle
c$$$              phaseatom = exp( img* 2d0*pi*sum(dq*pos(:,ibas1)) )
c$$$              m2   = m_indx    (ia2)
c$$$              l2   = l_indx    (ia2)
c$$$              n2   = n_indx    (ia2) + nc_max(l2,ibas2)
c$$$              lm2= l2**2 +l2+1 + m2
c$$$cccccccccccccccccccccccccccc
c$$$c Norm check test.
c$$$c          do j1= j1min,j1max
c$$$c          do j2= j2min,j2max
c$$$c            if(ia1==ia2) uum(j1,j2,ispin) = uum(j1,j2,ispin)
c$$$c     &        + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2)
c$$$c          enddo
c$$$c          enddo
c$$$ccccccccccccccccccccccccccccc
c$$$              do lm3= (l1-l2)**2+1, (l1+l2+1)**2 ! l3 can take |l1-l2|,...l1+l2
c$$$                l3 = ll(lm3)
c$$$ccccccccccccccccccccccccccccccc
c$$$c          ylk=0d0;   if(l3==0) ylk=1d0/sqrt(4*pi) !Y_00 only test
c$$$ccccccccccccccccccccccccccccccc
c$$$                ylk= cy(lm3)*yl(lm3)
c$$$                ppj = ppbrd(l1,n1,l2,n2,l3,ispin,ibas1) *cg(lm1,lm2, lm3)
c$$$     &          * fpi* img**l3* phaseatom * ylk
c$$$
c$$$! cg(lm1,lm2,lm3)= \int Y_lm3(\hat(r)) Y_lm2(\hat(r)) Y_lm1(\hat(r)) \frac{d \Omega}{4\pi}
c$$$! This is based on inverse expansion. See Rose.Eq.3.8.
c$$$                do j1= j1min,j1max
c$$$                  do j2= j2min,j2max
c$$$                    uum(j1,j2,ispin) = 
c$$$     &      uum(j1,j2,ispin) + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2) * ppj
c$$$                  enddo
c$$$                enddo
c$$$              enddo
c$$$ 1010       continue
c$$$ 1020     continue
c$$$
c$$$c ... Interstitial Plane Wave part
c$$$          call readgeig(q1, ngpmx, ispin, quu, geig1)
c$$$          call checkagree(q1,quu,' q1 ne quu eig')
c$$$          call readgeig(q2, ngpmx, ispin, quu, geig2)
c$$$          call checkagree(q2,quu,' q1 ne quu eig')
c$$$          do j1=j1min,j1max
c$$$            do j2=j2min,j2max
c$$$              uum(j1,j2,ispin)= uum(j1,j2,ispin) +
c$$$     &    sum( dconjg(geig1(1:ngp1,j1)) * matmul(ppovl,geig2(1:ngp2,j2)) )
c$$$            enddo
c$$$          enddo
c$$$        deallocate(cphi1,cphi2)
c$$$ 1050   continue
c$$$
c$$$c--- write resutlt
c$$$        write(6,*)' ============ result --- diagonal --- =============='
c$$$        do ispin = 1,nspin
c$$$          do j1=j1min,j1max
c$$$            do j2=j2min,j2max
c$$$              if(j1==j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)") 
c$$$     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
c$$$            enddo
c$$$          enddo
c$$$          write(6,*)'--- off diagonal ----------------------'
c$$$          do j1=j1min,j1max
c$$$            do j2=j2min,j2max
c$$$              if(j1/=j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)") 
c$$$     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$        write(6,*) ' ====== end ========================================'
c$$$c      stop ' ====== end ========================================'
c$$$        end
c$$$
c$$$        subroutine checkagree(a,b,char)
c$$$        real(8):: a(3),b(3)
c$$$        character*(*) :: char
c$$$        if(sum(abs(a-b))>1d-6) then
c$$$          write(6,*)' Error in checkagree:',char
c$$$Cstop2rx 2013.08.09 kino          stop ' Error in checkagree:'
c$$$          call rx( ' Error in checkagree:')
c$$$        endif
c$$$        end
