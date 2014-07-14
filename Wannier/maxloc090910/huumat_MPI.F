      program h_uumatrix 

c-------------------------------------------------
c  Calculate <u|u> matrix . u_kj(r) is the perodic part of eigencuntion.
c
c ixc=1: <u(k) | u(q_in)>
c ixc=2: <u(k) | u(k+b)>
c ixc=3: <u(k) | u(k+q0)>
c ixc=4: <u(k) | u(k+b)> for 1D (k=unshifted, shidted, or unshifted+q0)
c ixc=5: <u(q) | u(k)> (q from SYML, k on the FBZ mesh, closest to q)
c ixc=6: <u(q) | u(k)> (q from SYML, k on the FBZ mesh, 4 closest points)
c
c Takashi Miyake, Mar 2008, parallelized
c modified by Takahi Miyake, Aug, 2005 (ixc=4)
c modified by Takahi Miyake, Jul, 2005 (ixc=3)
c modified by Takahi Miyake, May, 2004
c originally written by Takao Kotani, April, 2004
c-------------------------------------------------
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
c RS: MPI module
      use rsmpi

      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     & , ngvecpf1(:,:), ngvecpf2(:,:), 
     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:),
     &phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:),
     &qbzs(:,:),qbz2(:,:)
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
     &, iopen,ico,lxd,lx !,ll
      real(8):: ylk
#ifdef COMMONLL
      integer(4):: ll(51**2)
      common/llblock/ll
#else
      integer(4) ll
#endif  
ctm
      integer(4) :: ixc,idummy,idummy2,i1,i2,i3,nbbloop,
     &              ifq0p,ifuu(2),
     &              ifbb,nbb,iko_ixs(2),iko_fxs(2),noxs(2),
     &              iqibz,iqbz,ibb,itmp,itmp2,iti,itf,
     &              nqibz2,nqbz2,iqb,ibb2,iqtmp,ibbtmp,ndg(3),
     &              nb1d,iq0i,nq
      integer(4),allocatable:: ikidx(:),ikbidx(:,:),
     &              ibidx(:,:),ibidxs(:,:),ibidx0(:,:,:),
     &              ij1idx(:),ij2idx(:)
      real(8),allocatable :: bbv(:,:)
      real(8) :: uunorm,dqx(3),dqx0(3),dq0(3),dg(3),dqmin(3),adq0,
     &           q0wf(3),wgt
      integer(4) :: iqindx
ctm, ixc=5
      integer(4) :: nqbandx,nqband,j1min_c(2),j1max_c(2),nbmin,nbmax,
     &        nmin,nmax,iq2,ntmp
      real(8) :: emin,emax,rydberg
      real(8),allocatable :: qbandx(:,:),qband(:,:),eband(:)
      logical:: lbnds
ctm, ixc=6
      integer(4),allocatable:: iqk4(:,:)
      real(8),allocatable:: wqk4(:,:)

ctm, MPI
!      include 'mpif.h'
!      integer(4):: ierr,istatus(MPI_STATUS_SIZE),iclose,input3(3)
      integer(4):: ierr,iclose,input3(3)
      integer(4):: myproc,nproc,nproc1,nproc2,nq_proc,ii,jj,kk,iftmp
      integer(4),allocatable:: iq_proc(:)
      character*8 xt
      character(4) charnum4
c-------------------------

ctm, MPI
c initialize MPI
!      call mpi_init(ierr)
!      call mpi_COMM_RANK(MPI_COMM_WORLD, myproc, ierr)
!      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

c RS: program begins here
      call RSMPI_Init()
      nproc = nproc_rsmpi
      myproc = myrank_rsmpi

ctm
c mode switch. --------------
      if (myproc.eq.0) then
      write(6,*) ' --- Choose omodes below -------------------'
      write(6,*) '  (1) (q,Qin),  (2) (q,q+b), (3) (q,q+q0)' 
      write(6,*) '  or (4) (q,q+b) for 1D' 
      write(6,*) ' --- Put number above ! --------------------'
      call readin5(ixc,idummy,idummy2)
      input3(1) = ixc
      input3(2) = idummy
      input3(3) = idummy2
      write(6,*) ' ixc=',ixc
      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
      endif ! myproc

      call MPI_Bcast(input3,3,MPI_INTEGER,0,
     &     MPI_COMM_WORLD,ierr)
      ixc=input3(1)
      idummy=input3(2)
      idummy2=input3(3)

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
      allocate(qbzs(3,nqbz))
      if(bzcase()==2.and.qbzreg()) then
      deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
      do i=1,nqbz
        qbzs(:,i) = qbz(:,i) -deltaq/2d0
        if (myproc.eq.0) write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
      enddo
ctm 050802
      if (ixc.le.6) stop "huumat: bzcase=2 is not implemented."
      endif
      if (myproc.eq.0) then
      print *,' ======== nqbz qbz  =',nqbz
      print *,' ======== nqibz ngrp=',nqibz,ngrp
      print *,  qbz 
      print *,'============================'
      print *
      endif

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

! read EF from 'BNDS' if exists
      if (ixc.eq.5.or.ixc.eq.6) then
      lbnds=.false.
      inquire(file='BNDS',exist=lbnds)      
      if (lbnds) then
        write(*,*)'Read EF from BNDS'
        open(99,file='BNDS',status='old')
        read(99,*)ntmp,ef
        close(99)
      else ! lbnds 
        call RSMPI_Stop( ' huumat_MPI: Cannot find BNDS')
      endif ! lbnds  
      endif ! ixc   

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
      if (myproc.eq.0) print *,' ngpmx=',ngpmx

c --- read radial functions PHIVC   (taken from hasfp0)
      if (myproc.eq.0) print *,' Go to readining phivc'
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
      if (myproc.eq.0) print *,' end of allocation'
      read(ifphi) nrad(1:nbas)
      read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
      nc_max=0
      do ibas=1,nbas
        ic = ibas
        if (myproc.eq.0) 
     &  write(6,*)' --- read PHIVC of ibas nrad=',ibas,nrad(ic)
        read(ifphi) ncore(ic), ncoremx                            !core
        read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
        if (myproc.eq.0) print *,' xxx0'
        read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)

        if (myproc.eq.0) 
     &  print *, 'ic icx=',ic,icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
        if(ic/=icx) then
           stop ' h_uu: ic/=icx'
        endif
        if (myproc.eq.0) 
     &  print *,' xxx1 ncoremx ncore(ic)=',ncoremx,ncore(ic)
        read(ifphi) rr(1:nrofi(ic),ic)
        if (myproc.eq.0) 
     &  print *,' xxx2 ncoremx ncore(ic)=',ncoremx,ncore(ic)

        if (myproc.eq.0) print *,' xxx2 nspin=',nspin
        rmax(ic) = rr(nrofi(ic),ic)
        do isp = 1, nspin
          if (myproc.eq.0) 
     & write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
          do ico = 1, ncore(ic) !core
            l =  lcindx(ico,ic)
            n =  ncindx(ico,ic)
            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
            phitotr(1:nrofi(ic),l,n, ic,isp)=              !core raw= core orthgonal
     &      phitoto(1:nrofi(ic),l,n, ic,isp)               !
            if(n>nc_max(l,ic)) nc_max(l,ic)=n
            if (myproc.eq.0) print *,' sss1c=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
          enddo
          do irad = 1, nrad(ic)   !valence
            l = lindx_r (irad,ic)
            n = nindx_r (irad,ic) + nc_max(l,ic)
            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal 
            read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
            if (myproc.eq.0) then
            print *,' sss1=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
            print *,' sss2=',sum(abs(phitotr(1:nrofi(ic),l,n, ic,isp)))
            endif
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
      allocate( geig1(ngpmx,nband),geig2(ngpmx,nband))
      if (myproc.eq.0) print *, 'end of initialization'

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
      if (myproc.eq.0) print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
!      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
      if (myproc.eq.0) write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
      nq0it = nq0i
      allocate( q0i(1:3,1:nq0i) ) !wqt(1:nq0i),
!      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      nq0ix = nq0i
      do i=1,nq0i
      read (101,* ) xxx,q0i(1:3,i)
      if(xxx==0d0 ) nq0ix = i-1
      enddo
      nq0i = nq0ix ! New nq0i July 2001
      if (myproc.eq.0) then
      write(6,*) ' Used k number in Q0P =', nq0i
      write(6,"(i3,2x, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
      endif
      close(101)

ctm 
c read BBVEC
      ifbb = iopen('BBVEC',1,0,0)
      if (ixc.le.3.or.ixc.eq.5.or.ixc.eq.6) then
        read(ifbb,*)
        read(ifbb,*)nbb,nqbz2
        if (nqbz.ne.nqbz2) stop 'readbb: nqbz is wrong!'
        allocate (bbv(3,nbb),ikbidx(nbb,nqbz))
        call  readbb(ifbb,nqbz,nspin,nbb,
     o             bbv,
     o             ikbidx,
     o             iko_ixs,iko_fxs,noxs)
      elseif (ixc.eq.4) then
        read(ifbb,*)
        read(ifbb,*)nbb,nqbz2,nb1d
        if (nbb.ne.1) stop 'readbb: nbb is wrong!'
        if (nqbz.ne.nqbz2) stop 'readbb: nqbz is wrong!'
        allocate (bbv(3,nbb),qbz2(3,nqbz))
        allocate(ibidx(nb1d,nqbz),ibidxs(nb1d,nqbz),
     >           ibidx0(nb1d,nqbz,nq0i),ij1idx(nb1d),ij2idx(nb1d))
        call  readbb1D(bzcase(),ifbb,nqbz,nspin,nbb,nb1d,nq0i,
     i             qbz,qbzs,q0i,ginv,
     o             bbv,
     o             ibidx,ibidxs,ibidx0,
     o             iko_ixs,iko_fxs,noxs)
      endif  

c open output files
c mode 1
      if (ixc .eq. 1) then

      ifq0p = iopen('QPwf',1,0,0)      
      read(ifq0p,*)
      read(ifq0p,"(d24.16,3x, 3d24.16)" ) (wgt,q0wf(1:3))

      if (myproc.eq.0) then
      ifuu(1) = iopen('UU0U.'//charnum4(0),0,-1,0)
      write(ifuu(1))'nqbz,nbb,iko_ixs(1),iko_fxs(1)'
      write(ifuu(1))nqbz,nbb,iko_ixs(1),iko_fxs(1)
      ifuu(1) = iclose('UU0U.'//charnum4(0))
      if (nspin.eq.2) then
        ifuu(2) = iopen('UU0D.'//charnum4(0),0,-1,0)
        write(ifuu(2))'nqbz,nbb,iko_ixs(2),iko_fxs(2)'
        write(ifuu(2))nqbz,nbb,iko_ixs(2),iko_fxs(2)
        ifuu(2) = iclose('UU0D.'//charnum4(0))
      endif
      endif

c mode 2
      elseif (ixc .eq. 2) then

      if (myproc.eq.0) then
      ifuu(1) = iopen('UUU.'//charnum4(0),0,-1,0)
      write(ifuu(1))'nqbz,nbb,iko_ixs(1),iko_fxs(1)'
      write(ifuu(1))nqbz,nbb,iko_ixs(1),iko_fxs(1)
      ifuu(1) = iclose('UUU.'//charnum4(0))
      if (nspin.eq.2) then
        ifuu(2) = iopen('UUD.'//charnum4(0),0,-1,0)
        write(ifuu(2))'nqbz,nbb,iko_ixs(2),iko_fxs(2)'
        write(ifuu(2))nqbz,nbb,iko_ixs(2),iko_fxs(2)
        ifuu(2) = iclose('UUD.'//charnum4(0))
      endif
      endif

c mode 3
      elseif (ixc .eq. 3) then

      if (myproc.eq.0) then
      ifuu(1) = iopen('UUq0U.'//charnum4(0),0,-1,0)
      write(ifuu(1))'nqbz,nq0i,iko_ixs(1),iko_fxs(1)'
      write(ifuu(1))nqbz,nq0i,iko_ixs(1),iko_fxs(1)
      ifuu(1) = iclose('UUq0U.'//charnum4(0))
      if (nspin.eq.2) then
        ifuu(2) = iopen('UUq0D.'//charnum4(0),0,-1,0)
        write(ifuu(2))'nqbz,nq0i,iko_ixs(2),iko_fxs(2)'
        write(ifuu(2))nqbz,nq0i,iko_ixs(2),iko_fxs(2)
        ifuu(2) = iclose('UUq0D.'//charnum4(0))
      endif
      endif

c mode 4
      elseif (ixc .eq. 4) then
         nbbloop = nq0i+1
         if (bzcase()==2) nbbloop = nq0i+2
      if (myproc.eq.0) then
      ifuu(1) = iopen('UU1dU.'//charnum4(0),0,-1,0)
      write(ifuu(1))nqbz,nq0i,nb1d,nbbloop
      ifuu(1) = iclose('UU1dU.'//charnum4(0))
      if (nspin.eq.2) then
        ifuu(2) = iopen('UU1dD.'//charnum4(0),0,-1,0)
        write(ifuu(2))nqbz,nq0i,nb1d,nbbloop
        ifuu(2) = iclose('UU1dD.'//charnum4(0))
      endif
      endif

c mode 5
      elseif (ixc .eq. 5) then

        nqbandx = 1000
        allocate(qbandx(3,nqbandx))
        call read_syml(qbandx,nqbandx,nqband)
        allocate(qband(3,nqband))
        qband(1:3,1:nqband) = qbandx(1:3,1:nqband)
        deallocate(qbandx)

      if (myproc.eq.0) then
      ifuu(1) = iopen('UUqkU.'//charnum4(0),0,-1,0)
      write(ifuu(1))'nqband,nqbz,nbb,iko_ixs(1),iko_fxs(1),nband'
      write(ifuu(1))nqband,nqbz,1,iko_ixs(1),iko_fxs(1),nband
      ifuu(1) = iclose('UUqkU.'//charnum4(0))
      if (nspin.eq.2) then
        ifuu(2) = iopen('UUqkD.'//charnum4(0),0,-1,0)
        write(ifuu(2))'nqband,nqbz,nbb,iko_ixs(2),iko_fxs(2),nband'
        write(ifuu(2))nqband,nqbz,1,iko_ixs(2),iko_fxs(2),nband
        ifuu(2) = iclose('UUqkD.'//charnum4(0))
      endif
      endif

c mode 6
      elseif (ixc .eq. 6) then

        nqbandx = 1000
        allocate(qbandx(3,nqbandx))
        call read_syml(qbandx,nqbandx,nqband)
        allocate(qband(3,nqband))
        qband(1:3,1:nqband) = qbandx(1:3,1:nqband)
        deallocate(qbandx)

        allocate(wqk4(4,nqband),iqk4(4,nqband))
        call findk4(qband,qbz,plat,qbas,nqband,nqbz,wqk4,iqk4)
!        do i = 1,nqband
!          write(*,"('***',i5,3f12.6)")i,qband(1:3,i)
!          do j = 1,4
!            write(*,"(i5,4f12.6)")iqk4(j,i),qbz(1:3,iqk4(j,i)),wqk4(j,i)
!          enddo   
!        enddo   
!        call RSMPI_Stop( ' huumat_MPI: debug')

      if (myproc.eq.0) then
      ifuu(1) = iopen('UUqk4U.'//charnum4(0),0,-1,0)
      write(ifuu(1))'nqband,nqbz,nbb,iko_ixs(1),iko_fxs(1),nband'
      write(ifuu(1))nqband,nqbz,4,iko_ixs(1),iko_fxs(1),nband
      ifuu(1) = iclose('UUqk4U.'//charnum4(0))
      if (nspin.eq.2) then
        ifuu(2) = iopen('UUqk4D.'//charnum4(0),0,-1,0)
        write(ifuu(2))'nqband,nqbz,nbb,iko_ixs(2),iko_fxs(2),nband'
        write(ifuu(2))nqband,nqbz,4,iko_ixs(2),iko_fxs(2),nband
        ifuu(2) = iclose('UUqk4D.'//charnum4(0))
      endif
      endif
c end of ixc-if
      endif


c======================================================================
c --- Set q1(j1range) q2(j2range)
c======================================================================
! Note that the true q when we generate eigenfunctions are q1x and q2x.
! q1-q1x should be a G vector.
! So you may need to take into account the phase shift to <u|u> vectors.
!
! --- I inserted checkagree to make sure that q1=q1x and q2=q2x ... 
!
ctm      q1 = qbz(:, 12)           ;  j1min=1;   j1max=8
ctm      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
c      q1 = qbz(:, 5);  j1min=1;   j1max=8
c      q2 = qbz(:, 5);  j2min=1;   j2max=8
c      q1 = qbz(:, 12)+q0i(1:3,1);  j1min=1;   j1max=8
c      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
ctm
      j1min = iko_ixs(1)
      j1max = iko_fxs(1)
      if (nspin .eq. 2) then
         if (iko_ixs(2) .lt. j1min) j1min = iko_ixs(2)
         if (iko_fxs(2) .gt. j1max) j1max = iko_fxs(2)
      endif
      j2min = j1min
      j2max = j1max
      if (ixc.eq.5.or.ixc.eq.6) then
         j1min = 1
         j1max = nband
      endif 
c======================================================================

      allocate( uum(j1min:j1max,j2min:j2max,nspin) )

ctm, MPI
      if (ixc.le.4) then
        nq = nqbz
      else ! ixc==5.or.ixc==6
        nq = nqband
      endif ! ixc
      nq_proc = nq / nproc
      nproc1 = nq - nq_proc*nproc
      nproc2 = nproc - nproc1
      allocate(iq_proc(nq_proc+1))
      iq_proc = 0
      iftmp = 100 + myproc
      write(iftmp,*)'*** myproc,i,q(i)'
      kk = 0
      do jj = 1,nproc
      do ii = 1,nq_proc+1
        if (jj.gt.nproc1.and.ii.gt.nq_proc) cycle
        kk = kk + 1
        if (myproc.eq.jj-1) then
          iq_proc(ii) = kk
           write(iftmp,*)myproc,ii,iq_proc(ii)
        endif  
      enddo ! ii
      enddo ! jj

      do 1070 ii = 1,nq_proc+1
         iqbz = iq_proc(ii)
         if (iqbz.eq.0) cycle
         write(*,*)'iq =',iqbz, 'out of',nq
         write(iftmp,*)'iq =',iqbz, 'out of',nq
c fileopen
         if (ixc.eq.1) then
           ifuu(1) = iopen('UU0U.'//charnum4(iqbz),0,-1,0)
           if (nspin.eq.2) ifuu(2) = iopen('UU0D.'//charnum4(iqbz),0,-1,0)
         elseif (ixc.eq.2) then   
           ifuu(1) = iopen('UUU.'//charnum4(iqbz),0,-1,0)
           if (nspin.eq.2) ifuu(2) = iopen('UUD.'//charnum4(iqbz),0,-1,0)
         elseif (ixc.eq.3) then   
           ifuu(1) = iopen('UUq0U.'//charnum4(iqbz),0,-1,0)
           if (nspin.eq.2) ifuu(2) = iopen('UUq0D.'//charnum4(iqbz),0,-1,0)
         elseif (ixc.eq.4) then   
           ifuu(1) = iopen('UU1dU.'//charnum4(iqbz),0,-1,0)
           if (nspin.eq.2) ifuu(2) = iopen('UU1dD.'//charnum4(iqbz),0,-1,0)
         elseif (ixc.eq.5) then   
           ifuu(1) = iopen('UUqkU.'//charnum4(iqbz),0,-1,0)
           if (nspin.eq.2) ifuu(2) = iopen('UUqkD.'//charnum4(iqbz),0,-1,0)
         elseif (ixc.eq.6) then   
           ifuu(1) = iopen('UUqk4U.'//charnum4(iqbz),0,-1,0)
           if (nspin.eq.2) ifuu(2) = iopen('UUqk4D.'//charnum4(iqbz),0,-1,0)
         endif   
ctm
c      do 1070 iqbz = 1,nq
c         write(*,*)'iq =',iqbz, 'out of',nq

      if (ixc.eq.1) then
         nbbloop = 1
      elseif (ixc.eq.2) then
         nbbloop = nbb
      elseif (ixc.eq.3) then
         nbbloop = nq0i
      elseif (ixc.eq.4) then
         nbbloop = nq0i+1
         if (bzcase()==2) nbbloop = nq0i+2
      elseif (ixc.eq.5) then
         nbbloop = 1
      elseif (ixc.eq.6) then
         nbbloop = 4
      endif

      do 1080 ibb = 1,nbbloop

         if (ixc.eq.1) then
            q1(:) = qbz(:,iqbz)
            q2(:) = q0wf(:)

         elseif (ixc.eq.2) then

            if (myproc.eq.0)write(*,*)'ibb  =',ibb, 'out of',nbb

            iqb = ikbidx(ibb,iqbz)
            q1(:) = qbz(:,iqbz)
c            q2 = qbz(:,iqb)
            q2(:) = q1(:) + bbv(:,ibb)

            if (iqb .lt. iqbz) then
               iqtmp = iqb
               do ibb2 = 1,nbb
                  itmp = ikbidx(ibb2,iqtmp)
                  if (itmp .eq. iqbz) then
                     ibbtmp = ibb2
                     goto 1200
                  endif
               enddo
               stop 'huumat: (iq,ib) error'
 1200          continue

               do ispin = 1,nspin
                  write(ifuu(ispin))-20
                  write(ifuu(ispin))iqbz,ibb,iqtmp,ibbtmp
               enddo
               goto 1090
            endif

         elseif (ixc.eq.3) then

            if (myproc.eq.0)write(*,*)'ibb  =',ibb, 'out of',nbb

            q1(:) = qbz(:,iqbz)
c            q2 = qbz(:,iqb)
            q2(:) = q1(:) + q0i(:,ibb)

         elseif (ixc.eq.4) then

            if (myproc.eq.0)write(*,*)'iq0i  =',ibb, 'out of',nbbloop

            if (bzcase()==2) then 
              if (ibb.eq.1) then
                q1(:) = qbz(:,iqbz)
                q2(:) = q1(:) + bbv(:,1)
                ij1idx(:) = ibidx(:,iqbz)
                iqtmp = iqindx(q2(:), ginv,qbz,nqbz)
                ij2idx(:) = ibidx(:,iqtmp)
              elseif (ibb.eq.2) then  
                q1(:) = qbzs(:,iqbz)
                q2(:) = q1(:) + bbv(:,1)
                ij1idx(:) = ibidxs(:,iqbz)
                iqtmp = iqindx(q2(:), ginv,qbzs,nqbz)
                ij2idx(:) = ibidxs(:,iqtmp)
              else
                iq0i = ibb-2
                q1(:) = qbz(:,iqbz)+q0i(:,iq0i)
                q2(:) = q1(:) + bbv(:,1)
                ij1idx(:) = ibidx0(:,iqbz,iq0i)
                do iq = 1,nqbz
                   qbz2(:,iq) = qbz(:,iq) + q0i(:,iq0i)
                enddo   
                iqtmp = iqindx(q2(:), ginv,qbz2,nqbz)
                ij2idx(:) = ibidx0(:,iqtmp,iq0i)
              endif   
            else  
              if (ibb.eq.1) then
                q1(:) = qbz(:,iqbz)
                q2(:) = q1(:) + q0i(:,ibb)
                ij1idx(:) = ibidx(:,iqbz)
                iqtmp = iqindx(q2(:), ginv,qbz,nqbz)
                ij2idx(:) = ibidx(:,iqtmp)
              else
                iq0i = ibb-1
                q1(:) = qbz(:,iqbz)+q0i(:,iq0i)
                q2(:) = q1(:) + bbv(:,1)
                ij1idx(:) = ibidx0(:,iqbz,iq0i)
                do iq = 1,nqbz
                   qbz2(:,iq) = qbz(:,iq) + q0i(:,iq0i)
                enddo   
                iqtmp = iqindx(q2(:), ginv,qbz2,nqbz)
                ij2idx(:) = ibidx0(:,iqtmp,iq0i)
              endif   
            endif  
         elseif (ixc.eq.5) then
            q1(:) = qband(:,iqbz)
            call findk(q1,qbz,plat,qbas,nqbz,iq2)
            q2 = qbz(:,iq2)
            if (myproc.eq.0)write(*,*)'q2  =',q2
            call getkeyvalue("GWinput","uu_emin",emin,default=-999d0 )
            call getkeyvalue("GWinput","uu_emax",emax,default=999d0 )
            call getkeyvalue("GWinput","uu_nbmin",nbmin,default=-999 )
            call getkeyvalue("GWinput","uu_nbmax",nbmax,default=999 )
            allocate(eband(nband))
            do is = 1,nspin
              j1min_c(is) = j1min
              j1max_c(is) = j1max
              call readeval(q1,is,eband)
              eband = (eband-ef)*rydberg()
              nmin = nband
              do i = nband,1,-1
               if (eband(i) .ge. emin) nmin = i
              enddo
              nmin = max0(nmin,nbmin)
              nmax = 1
              do i = 1,nband
               if (eband(i) .le. emax) nmax = i
              enddo
              nmax = min0(nmax,nbmax)
              if (nmin.gt.j1min_c(is)) j1min_c(is) = nmin
              if (nmax.lt.j1max_c(is)) j1max_c(is) = nmax
              if (j1min_c(is).gt.j1max_c(is))stop' huumat: window error'
              if (iqbz.eq.1.and.ibb.eq.1) then
               write(*,*)'uu_emin',emin
               write(*,*)'uu_emax',emax
               write(*,*)'uu_nbmin',nbmin
               write(*,*)'uu_nbmin',nbmax
               write(*,*)'j1min_c',j1min_c(is)
               write(*,*)'j1max_c',j1max_c(is)
               write(*,*)'nband  ',nband
              endif   
            enddo  
            deallocate(eband)
         elseif (ixc.eq.6) then
            q1(:) = qband(:,iqbz)
            iq2 = iqk4(ibb,iqbz)
            q2 = qbz(:,iq2)
            if (myproc.eq.0)write(*,*)'q2  =',q2
            call getkeyvalue("GWinput","uu_emin",emin,default=-999d0 )
            call getkeyvalue("GWinput","uu_emax",emax,default=999d0 )
            call getkeyvalue("GWinput","uu_nbmin",nbmin,default=-999 )
            call getkeyvalue("GWinput","uu_nbmax",nbmax,default=999 )
            allocate(eband(nband))
            do is = 1,nspin
              j1min_c(is) = j1min
              j1max_c(is) = j1max
              call readeval(q1,is,eband)
              eband = (eband-ef)*rydberg()
              nmin = nband
              do i = nband,1,-1
               if (eband(i) .ge. emin) nmin = i
              enddo
              nmin = max0(nmin,nbmin)
              nmax = 1
              do i = 1,nband
               if (eband(i) .le. emax) nmax = i
              enddo
              nmax = min0(nmax,nbmax)
              if (nmin.gt.j1min_c(is)) j1min_c(is) = nmin
              if (nmax.lt.j1max_c(is)) j1max_c(is) = nmax
              if (j1min_c(is).gt.j1max_c(is))stop' huumat: window error'
              if (iqbz.eq.1.and.ibb.eq.1) then
               write(*,*)'uu_emin',emin
               write(*,*)'uu_emax',emax
               write(*,*)'uu_nbmin',nbmin
               write(*,*)'uu_nbmin',nbmax
               write(*,*)'j1min_c',j1min_c(is)
               write(*,*)'j1max_c',j1max_c(is)
               write(*,*)'nband  ',nband
              endif   
            enddo  
            deallocate(eband)
c end of ixc-if
         endif

C --- q1x and q2x
      call readqg0('QGpsi',q1,ginv,  q1x, ngp1)
      call readqg0('QGpsi',q2,ginv,  q2x, ngp2)
c      call checkagree(q1,q1x,' q1 ne q1x')
c      call checkagree(q2,q2x,' q2 ne q2x')
!      write(6,"('q1 q1x=',3f9.4,3x,3f9.4)") q1,q1x
!      write(6,"('q2 q2x=',3f9.4,3x,3f9.4)") q2,q2x

c ... dq 
      dqx = q1x-q2x
      call q2q0(dqx,plat,dqx0)
      call q02q0g0(dqx0,dq0,ndg)
      ndg = -ndg
      call q02q(dq0,qbas,dq)
      dg = dq - dqx

      if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)

      absdq = sqrt(sum(dq**2))
      absqg2 = (2*pi/alat)**2 *sum(dq**2)
      absqg =sqrt(absqg2)

C --- ppovl= <P_{q1+G1}|P_{q2+G2}> 
      allocate( ngvecpf1(3,ngp1), ngvecpf2(3,ngp2), ppovl(ngp1,ngp2) )
      call readqg('QGpsi',q1,ginv, q1x, ngp1, ngvecpf1)
      call readqg('QGpsi',q2,ginv, q2x, ngp2, ngvecpf2)
c      call checkagree(q1,q1x,' q1 ne q1x xxx2')
c      call checkagree(q2,q2x,' q2 ne q2x xxx2')
c      print *,' ngp1,ngp2 sum check=',ngp1,ngp2,sum(abs(ngvecpf1)),sum(abs(ngvecpf2))
c 19June2004 T.Miyake  ccccccccccccccccccccccccccccccc
      do i = 1,3
      ngvecpf2(i,1:ngp2) = ngvecpf2(i,1:ngp2) + ndg(i)
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mkppovl2(alat,plat,qbas, !22April2004 
     &    ngp1, ngvecpf1, 
     &    ngp2, ngvecpf2, 
     &    nbas, rmax, pos,  
     o    ppovl)
c      print *,' end of mkppovl2'

c ... lxx and allocations
      lxx=2*(nl-1)
      allocate( ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nspin,nbas),
     &   rprodx(nrx,0:lxx),
     &   phij(0:lxx),psij(0:lxx),rphiphi(nrx))

c ... dq 
c      dq = q1x-q2x
c      if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
c
c      absdq = sqrt(sum(dq**2))
c      absqg2 = (2*pi/alat)**2 *sum(dq**2)
c      absqg =sqrt(absqg2)

c ... YL(dq)
      allocate(cy((lxx+1)**2),yl((lxx+1)**2))
      call sylmnc(cy,lxx)
      call sylm(dq/absdq,yl,lxx,r2s) !spherical factor Y(dq)

C --- radial integral  ppbrd = <phi phi j_l>
      ppbrd=0d0
      do 900 ibas = 1,nbas
        ic = ibas
c        write(6,"(' nindx=',10i3)") nindx(1:nl,ic)
c        write(6,*)' radial integral ibas=',ibas
        do ir =1,nrofi(ic)
          call bessl(absqg2*rr(ir,ibas)**2,lxx,phij,psij)
c  phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
          do lx = 0, lxx
            if(rr(ir,ibas)==0d0) then
              rprodx(ir,lx)=0d0
            else   
              rprodx(ir,lx) = rr(ir,ibas)* phij(lx)* (absqg*rr(ir,ibas))**lx  
            endif  
            ! = r \times j_l(|dq|r)  !bessel function
          enddo
ccccccccccccccccccccccc
c          write(1100,"(10d13.6)")rr(ir,ibas),rprodx(ir,0:lxx) !,phij(0:lxx)
ccccccccccccccccccccccc
        enddo
        do 125 isp = 1,nspin
         do 25 l1 = 0, nl-1      
         do 25 n1 = 1, nindx(l1+1,ic)
         do 25 l2 = 0, nl-1      
         do 25 n2 = 1, nindx(l2+1,ic)
          rphiphi(1)       = 0d0
          rphiphi(2:nrofi(ic)) = phitoto(2:nrofi(ic),l1,n1,ic,isp)
     &                          *phitoto(2:nrofi(ic),l2,n2,ic,isp)/rr(2:,ic) ! phi = u = r \phi
         do 25 lx = 0, 2*(nl-1)
          if(lx <abs(l1-l2) .or. l1+l2<lx) cycle
          call gintxx( rprodx(1,lx), rphiphi,aa(ic),bb(ic),nrofi(ic),
     &        ppbrd(l1, n1,l2, n2, lx, isp,ibas) )
c          if(l1==l2.and.n1==n2.and.lx==0) 
c         print *, ' ppbrd=',l1,n1,ppbrd(l1, n1,l2, n2, lx, isp,ibas)
   25    continue
 125    continue
 900  continue

C --- Calcuate <u{q1x j1} | u_{q2x j2}> 
c              = < exp(i(q1x-q2x)r) psi^*{q1x j1} psi_{q2x j2} > 
c ... MT part
cr   ldim2 = nlmto
cr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
cr   l_indx   (1:ldim2) : l index
cr   ibas_indx(1:ldim2) : ibas index.
      uum = 0d0
      do 1050 ispin=1,nspin
       allocate(cphi1 (nlmto,nband),cphi2(nlmto,nband) )
       call readcphi(q1, nlmto, ispin, quu, cphi1)
c       call checkagree(q1,q1x,' q1 ne quu')
       call readcphi(q2, nlmto, ispin, quu, cphi2)
c       call checkagree(q2,quu,' q2 ne quu')

       do 1020 ia1 = 1,nlmto
         ibas1= ibas_indx(ia1)
         l1   = l_indx    (ia1)
         m1   = m_indx    (ia1)
         n1   = n_indx    (ia1) + nc_max(l1,ibas1)
         lm1  = l1**2+l1+1  + m1
       do 1010 ia2 = 1,nlmto
         ibas2 = ibas_indx(ia2)
         if(ibas2/=ibas1) cycle
         phaseatom = exp( img* 2d0*pi*sum(dq*pos(:,ibas1)) )
         m2   = m_indx    (ia2)
         l2   = l_indx    (ia2)
         n2   = n_indx    (ia2) + nc_max(l2,ibas2)
         lm2= l2**2 +l2+1 + m2
cccccccccccccccccccccccccccc
c Norm check test. 
c          do j1= j1min,j1max
c          do j2= j2min,j2max
c            if(ia1==ia2) uum(j1,j2,ispin) = uum(j1,j2,ispin) 
c     &        + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2)
c          enddo
c          enddo
ccccccccccccccccccccccccccccc
         do lm3= (l1-l2)**2+1, (l1+l2+1)**2 ! l3 can take |l1-l2|,...l1+l2
          l3 = ll(lm3) 
ccccccccccccccccccccccccccccccc
c          ylk=0d0;   if(l3==0) ylk=1d0/sqrt(4*pi) !Y_00 only test
ccccccccccccccccccccccccccccccc
          ylk= cy(lm3)*yl(lm3)
          ppj = ppbrd(l1,n1,l2,n2,l3,ispin,ibas1) *cg(lm1,lm2, lm3)
     &          * fpi* img**l3* phaseatom * ylk

! cg(lm1,lm2,lm3)= \int Y_lm3(\hat(r)) Y_lm2(\hat(r)) Y_lm1(\hat(r)) \frac{d \Omega}{4\pi}
! This is based on inverse expansion. See Rose.Eq.3.8.
          if (ixc.ne.4) then
c            do j1= j1min,j1max
c            do j2= j2min,j2max
            do j1= iko_ixs(ispin),iko_fxs(ispin)
            do j2= iko_ixs(ispin),iko_fxs(ispin)
             uum(j1,j2,ispin) = 
     &       uum(j1,j2,ispin) + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2)
     &                        * ppj
            enddo
            enddo
          else
            do j1= iko_ixs(ispin),iko_fxs(ispin)
            do j2= iko_ixs(ispin),iko_fxs(ispin)
             uum(j1,j2,ispin) =  uum(j1,j2,ispin) 
     &  + dconjg(cphi1(ia1,ij1idx(j1)))*cphi2(ia2,ij2idx(j2)) * ppj
            enddo
            enddo
          endif   
         enddo
 1010  continue
 1020  continue

c ... Interstitial Plane Wave part
       call readgeig(q1, ngpmx, ispin, quu, geig1)
c       call checkagree(q1,quu,' q1 ne quu eig')
       call readgeig(q2, ngpmx, ispin, quu, geig2)
c       call checkagree(q2,quu,' q1 ne quu eig')
       if (ixc.eq.5.or.ixc.eq.6) then
         do j1= j1min_c(ispin),j1max_c(ispin)
         do j2= iko_ixs(ispin),iko_fxs(ispin)
          uum(j1,j2,ispin)= uum(j1,j2,ispin) +
     &      sum( dconjg(geig1(1:ngp1,j1))
     &    * matmul(ppovl,geig2(1:ngp2,j2)) )
         enddo
         enddo
       elseif (ixc.ne.4) then
c         do j1=j1min,j1max
c         do j2=j2min,j2max
         do j1= iko_ixs(ispin),iko_fxs(ispin)
         do j2= iko_ixs(ispin),iko_fxs(ispin)
          uum(j1,j2,ispin)= uum(j1,j2,ispin) +
     &      sum( dconjg(geig1(1:ngp1,j1))
     &    * matmul(ppovl,geig2(1:ngp2,j2)) )
         enddo
         enddo
       else
         do j1= iko_ixs(ispin),iko_fxs(ispin)
         do j2= iko_ixs(ispin),iko_fxs(ispin)
          uum(j1,j2,ispin)= uum(j1,j2,ispin) +
     &      sum( dconjg(geig1(1:ngp1,ij1idx(j1)))
     &    * matmul(ppovl,geig2(1:ngp2,ij2idx(j2))) )
         enddo
         enddo
       endif   
       deallocate(cphi1, cphi2)
 1050 continue

ctm
      do ispin = 1,nspin
         iti = iko_ixs(ispin)
         itf = iko_fxs(ispin)
         write(ifuu(ispin))-10
       if (ixc.eq.5) then
            write(ifuu(ispin)) iqbz,iq2
            write(ifuu(ispin))
     &         ((uum(j1,j2,ispin),j1=1,nband),j2=iti,itf)
       elseif (ixc.eq.6) then
            write(ifuu(ispin)) iqbz,ibb,iq2,wqk4(ibb,iqbz)
            write(ifuu(ispin))
     &         ((uum(j1,j2,ispin),j1=1,nband),j2=iti,itf)
       else ! ixc
         if (ixc.eq.1) then
            write(ifuu(ispin)) iqbz
         elseif (ixc.eq.2) then
            write(ifuu(ispin))
     &            iqbz,ibb,ikbidx(ibb,iqbz)
         elseif (ixc.eq.3.or.ixc.eq.4) then
            write(ifuu(ispin)) iqbz,ibb
         else
         endif
         write(ifuu(ispin))
     &         ((uum(j1,j2,ispin),j1=iti,itf),j2=iti,itf)
       endif ! ixc
      enddo

check write
!      if (mod(iqbz,10).eq.1) then
!         do ispin = 1,nspin
!            iti = iko_ixs(ispin)
!            itf = iko_fxs(ispin)
!            do j1 = iti,itf
!               uunorm = 0d0
!               do j2 = iti,itf
!                  uunorm = uunorm + abs(uum(j1,j2,ispin))**2
!               enddo
!               write(*,980)ispin,iqbz,ibb,j1,uunorm
!            enddo
!         enddo
!      endif

c--- write resutlt 
c      print *,' ============ result --- diagonal --- =============='
c      do ispin = 1,nspin
c      do j1=j1min,j1max
c      do j2=j2min,j2max
c       if(j1==j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)") 
c     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
c      enddo
c      enddo
c      print *,'--- off diagonal ----------------------'
c      do j1=j1min,j1max
c      do j2=j2min,j2max
c       if(j1/=j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)") 
c     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
c      enddo
c      enddo
c      enddo

      deallocate(ngvecpf1, ngvecpf2, ppovl, ppbrd,
     &   rprodx, phij, psij, rphiphi, 
     &   cy, yl)

 1090 continue

c end of bb-loop
 1080 continue

c fileclose
         if (ixc.eq.1) then
           ifuu(1) = iclose('UU0U.'//charnum4(iqbz))
           if (nspin.eq.2) ifuu(2) = iclose('UU0D.'//charnum4(iqbz))
         elseif (ixc.eq.2) then
           ifuu(1) = iclose('UUU.'//charnum4(iqbz))
           if (nspin.eq.2) ifuu(2) = iclose('UUD.'//charnum4(iqbz))
         elseif (ixc.eq.3) then
           ifuu(1) = iclose('UUq0U.'//charnum4(iqbz))
           if (nspin.eq.2) ifuu(2) = iclose('UUq0D.'//charnum4(iqbz))
         elseif (ixc.eq.4) then
           ifuu(1) = iclose('UU1dU.'//charnum4(iqbz))
           if (nspin.eq.2) ifuu(2) = iclose('UU1dD.'//charnum4(iqbz))
         endif

c end of iqbz-loop
 1070 continue


 980  format('uu: norm check',i4,i5,2i4,f12.6)
 990  format(3i6)

      if (myproc.eq.0)
     & print *, ' ====== end ========================================'
c      stop ' ====== end ========================================'
c finalize MPI
!      call mpi_finalize(ierr)
      call RSMPI_Finalize()
      end

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
      subroutine  readbb1D(ibzcase,ifbb,nqbz,nspin,nbb,nband,nq0i,
     i                   qbz,qbzs,q0i,ginv,
     o                   bbv,
     o                   ibidx,ibidxs,ibidx0,
     o                   iko_ixs,iko_fxs,noxs)

      implicit real*8(a-h,o-z)
      parameter (eps = 1d-4)

      real (8) :: u(3),bbv(3,nbb),ginv(3,3),qbz(3,nqbz),
     >            qbzs(3,nqbz),q0i(3,nq0i),q(3),qbz2(3,nqbz)
      integer (4) :: iopen,iqindx,itmp,i,
     i           iko_ixs(2),iko_fxs(2),noxs(2)
      integer(4):: ibidx(nband,nqbz),ibidxs(nband,nqbz)
     >            ,ibidx0(nband,nqbz,nq0i)
     >            ,ibidx2(nband)
c      integer(4),allocatable:: ikidx(:),ikbidx(:,:)

c read bb-vector
      do i = 1,nbb
         read(ifbb,*)bbv(1,i),bbv(2,i),bbv(3,i)
      enddo
c read q-vector

c unshifted grid
      do iq = 1,nqbz
        read(ifbb,*)q(1:3)
        read(ifbb,*)ibidx2(1:nband)
        itmp = iqindx(q(:), ginv,qbz,nqbz)
        ibidx(:,itmp) = ibidx2(:) 
      enddo   

c shifted grid
      if (ibzcase.eq.2) then
        do iq = 1,nqbz
          read(ifbb,*)q(1:3)
          read(ifbb,*)ibidx2(1:nband)
          itmp = iqindx(q(:), ginv,qbzs,nqbz)
          ibidxs(:,itmp) = ibidx2(:) 
        enddo   
      endif   

c +q0 grid
      do iq0=1,nq0i
        do iq=1,nqbz
          qbz2(:,iq) = qbz(:,iq)+q0i(:,iq0)
        enddo
        do iq = 1,nqbz
          read(ifbb,*)q(1:3)
          read(ifbb,*)ibidx2(1:nband)
          itmp = iqindx(q(:), ginv,qbz2,nqbz)
          ibidx0(:,itmp,iq0) = ibidx2(:) 
        enddo   
      enddo   

c iko_ixs and iko_fxs
      do is = 1,nspin
         iko_ixs(is) = 1
         iko_fxs(is) = nband
         noxs(is) = nband
      enddo

      return
      end
c-----------------------------------------------------------------------
