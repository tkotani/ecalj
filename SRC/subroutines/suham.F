      module m_suham
      integer, allocatable,protected ::  iv_a_oindxo (:)
      integer, allocatable,protected ::  iv_a_ooffH (:)
      integer, protected:: ham_npwmin, ham_ndham,pot_nlma, pot_nlml,ham_lmxax,ham_ldham(16)

      contains
      subroutine m_suham_init(sspec,ssite)
      use m_struc_def  
      use m_lmfinit,only:lso, ctrl_nbas,ctrl_nspec, ctrl_nl, ctrl_nspin,
     &     ham_pwmode,ham_pwemax,ham_pwemin,lat_alat, lat_alat,lat_tolft,nkaph
      use m_supot,only: lat_ng, rv_a_ogv
      use m_lattic,only: lat_qlat,lat_plat
C- Hamiltonian setup
C ----------------------------------------------------------------------
Cio sspec,ssite can be modified.
Cio  sspec 
Cio    Elts read: lmxa lmxl p pz
Cio    Stored:    *
Cio    Passed to: atfold makidx showbs sugcut
Cio  ssite :struct for site-specific information; see routine usite
Cio    Elts read: spec
Cio    Stored:    pnu pz
Cio    Passed to: showbs pvioeu
Cl   ndim  :total number of lmto orbitals = nl**2 * nbas
Cl   npwmin:lower limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   npwmax:upper limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   nqdiv: loop over mesh of q-points to estimate npwmax
Cl        : nqdiv is fineness of q-mesh.
Clxxx  npwpad: a 'safety' padding to npwmax in case npwmax
Clxxx        : underestimates actual upper limit
Cr Remarks
Cr   This routine generates energy-independent hamiltonian setup.
Cr  *It generates and packs a table of hamiltonian offsets offH,
Cr   orbital permutation indices oindxo.
Cr
Cu Updates
Cu   07 Jul 08 Make sham->ndham = estimate for upper dimension of
Cu             hamiltonian, including possible PW part
Cu             Make sham->lmxax = largest lmax in basis
Cu   18 Apr 05 Force small parameter p -> 0 in 2-center turned on (ASA)
Cu   14 Feb 03 Makes and packs magnetic field B.sigma
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
c      type(s_ctrl)::sctrl
c      type(s_lat)::slat
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)
c      type(s_ham)::sham
c      type(s_pot)::spot
c      type(s_bz)::sbz
c      character*(*) sstrn
C ... Local parameters
      integer nkap0,n0H
      parameter (nkap0=3,n0H=5)
      integer hord,i,i1,i2,iprint,ldham(16),ldim, !lgen3,lham,
     .lgunit,lidim,lihdim,nbas,nclasp,ndim,neul,nl,
c lncol,
     .nlspcp,nsp,nspc,nspec,nspx,nttab,partok,igets,nvi,nvl, !lfp,
     .ib,is,lmxa,lmxl,stdo,isw,nbf,lmxax,pwmode,
     .j1,j2,j3,m,npw,npwmin,npwmax,ndham
      double precision pwemin,pwemax,plat(3,3),qlat(3,3),q(3),Gmin,Gmax,
     .xx
      integer n0,nqdiv,npwpad
      parameter (n0=10,nqdiv=12)
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))

C ... Needed for iostr
      double precision ckbas,cksumf,kap2(10)
      integer:: oalph , oiax , ontab , os , ng
      real(8),allocatable :: rv_a_ogvx(:)

      integer ,allocatable :: kv_iv(:)
      integer ,allocatable :: igv2_iv(:)

      integer:: obas , oo , opp , oeula , opti , obs , omagf,i_copy_size
      logical ltmp,iostr,bittst,adnf
      character*80 outs,gfopts
      double precision qss(4),vmtz,pnu(n0,2),pz(n0,2),alat,tolgv
      integer,allocatable:: ips(:)

      nbas=ctrl_nbas
      nspec=ctrl_nspec
      nl=ctrl_nl
c      lncol=ctrl_lncol
      nsp=ctrl_nspin
c      nkaph = globalvariables%nkaph
      pwmode=ham_pwmode
      pwemin=ham_pwemin
      pwemax=ham_pwemax
c      npwpad=ham_npwpad

      nspc  = 1
      nspx  = nsp
c      if (lncol .ne. 0) then
      if (lso .ne. 0) then
        nspc = 2
        nspx = 1
      endif
      ndim = nbas * nl**2 * nkaph
      if (allocated(iv_a_ooffH)) deallocate(iv_a_ooffH)
      allocate(iv_a_ooffH(n0h*nkap0*(nbas+1)))
      iv_a_ooffH(:)=0
      allocate(iv_a_oindxo(abs(ndim+3)))
      stdo = lgunit(1)

C --- Hamiltonian offsets, orbital permutation table ---
      call iinit(ldham,16)
      if (mod(pwmode,10) .eq. 2) then
        call info0(20,0,0,' suham: LMTO basis will be excluded')
      else
        allocate(ips(nbas))
        do  ib = 1, nbas
          ips(ib) = ssite(ib)%spec
        enddo  
        call makidx ( nl , nkaph , 1 , nbas , 0 , sspec , ips !sarray%iv_a_oips
     .     , iv_a_ooffH , iv_a_oindxo , ldham )
        deallocate(ips)
      endif
C     Default value for hamiltonian dimension
      ham_ndham=ldham(1)
C ... Make rest of ldham
      ldham(4) = nspc
      ldham(5) = ldham(1) * nspc
      ldham(6) = ldham(2) * nspc
      ldham(7) = ldham(3) * nspc
      ldham(8) = nspx
      ham_ldham=ldham

C ------- Potential- and implementation- specific setup -------
      alat=lat_alat
      plat=lat_plat
      qlat=lat_qlat
      nvi = 0
      nvl = 0
      lmxax = -1
      do  20  ib = 1, nbas

        is=ssite(ib)%spec
        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl
        i_copy_size=size(sspec(is)%p)
        call dcopy(i_copy_size,sspec(is)%p,1,pnu,1)
        i_copy_size=size(sspec(is)%pz)
        call dcopy(i_copy_size,sspec(is)%pz,1,pz,1)
c        pnu(:)= sspec(is)%p(:)
c        pz(:) = sspec(is)%pz(:)
cccccccccccccc
c        print *,'suham 111 pnu=',i,pnu(1:10,1:2)
cccccccccccccc
        call dcopy(n0,pz,1,pz(1,2),1)
cc        pz(:,2)=pz
C         Augmentation dimensioning parameters
        nvi = nvi + (lmxa+1)**2
        nvl = nvl + (lmxl+1)**2
C         Poke spec starting pnu to site
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,pz,1,ssite(ib)%pz,1)
cc        ssite(ib)%pnu=pnu
cc        ssite(ib)%pz=pz
C         Find largest lmxa
        lmxax = max(lmxax,lmxa)
   20 continue
      ham_lmxax=lmxax
      pot_nlma=nvi
      pot_nlml=nvl
      alat=lat_alat
      call info5(30,0,0,' suham :  %i augmentation'//
     .' channels, %i local potential channels  Maximum lmxa=%i', nvi,nvl,lmxax,0,0)
c      ng=lat_ng
c      tolgv=lat_tolft
c      call sugcut ( 1 , nspec , sspec , alat , ng , rv_a_ogv , tolgv  )
      call sugcut ( 1 , sspec) !modify sspec%ngcut

C   ... PW setup : estimate upper bound to number of G vectors
C       to set up upper bound to hamiltonian dimension
      pwmode = ham_pwmode
      pwemax = ham_pwemax
      if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
        Gmin = dsqrt(pwemin)
        Gmax = dsqrt(pwemax)
        if (mod(pwmode/10,10) .eq. 1) then
          call info0(70,1,0,' Estimate max size of PW basis from'//
     .    'combinations of recip. lattice vectors ...')
          npwmax = -1
          npwmin = 99999
          do  j1 = 0, nqdiv
            do  j2 = 0, nqdiv
              do  j3 = 0, nqdiv
                do   m = 1, 3
                  q(m) = (qlat(m,1)/nqdiv)*j1 +
     .            (qlat(m,2)/nqdiv)*j2 +
     .            (qlat(m,3)/nqdiv)*j3
                enddo
                call pshpr(iprint()-40)
ctakao nqdiv is large but fixed... So this is safer?
c              call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
c     .          xx,xx,xx)
c             write(6,"('GVLST2:goto q=',3f13.5)")q
c             call shorbz(q,q,qlat,plat)
                call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
     &          xx,xx,xx)
                call poppr
                npwmin = min(npwmin,npw)
                npwmax = max(npwmax,npw)
              enddo
            enddo
          enddo
c          if (npwpad .lt. 0) then
            npwpad = max(nint((npwmax-npwmin)*0.2d0),3)
c          endif
        else
          call dpzero(q,3)
          call pshpr(iprint()-40)
          call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
     .    xx,xx,xx)
          call poppr
          npwmin = npw
          npwmax = npw
          npwpad = 0
        endif
        ndham = npwmax + npwpad
        if (mod(pwmode,10) .ne. 2) ndham = ldham(1) + npwmax + npwpad

        ham_npwmin=npwmin
c        ham_npwpad=npwpad
        ham_ndham=ndham

        if (mod(pwmode/10,10) .eq. 1) then
          call info2(20,1,0,' suham:  q-dependent PW basis with'//
     .    '  Emin = %d < E < %d.',pwemin,pwemax)
c          call info5(30,0,0,'%9fEst. min,max PW dimension = %i,%i.'//
c     .    '  Use npwpad = %i => ndham = %i',
c     .    npwmin,npwmax,npwpad,ndham,0)
        else
          call info5(20,0,0,' suham:  PW basis with  %d < E < '//
     .    '%d  =>  npw = %i,  ndham = %i',
     .    pwemin,pwemax,npw,ndham,0)
        endif

C    ...  Printout
        if (iprint() .ge. 40) then
          call pshpr(0)
          call dpzero(q,3)
          call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
     .    xx,xx,xx)
          call poppr
          call info0(40,1,-1,' G vectors at the Gamma point:')
          allocate(igv2_iv(3*npw))
          allocate(kv_iv(3*npw))
          allocate(rv_a_ogvx(abs(3*npw))); rv_a_ogvx(:)=0.0d0
          call pshpr(iprint())
          if (iprint() .ge. 50) call setpr(100)
          call gvlst2 ( alat , plat , q , 0 , 0 , 0 , gmin , gmax , 0 ,
     .     8 + 2 , npw , npw , kv_iv , rv_a_ogvx , xx , igv2_iv )
          call poppr
          if (allocated(kv_iv)) deallocate(kv_iv)
          if (allocated(igv2_iv)) deallocate(igv2_iv)
          deallocate(rv_a_ogvx)
        endif
      endif
      end subroutine m_suham_init
      end module m_suham
cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine gen_hamindex(siginit,ljobgw)
      use m_suham,only: iprmb=>iv_a_oindxo, ham_ldham
      use m_supot,only: rv_a_ogv
      use m_mksym,only: rv_a_osymgr,rv_a_oag,lat_nsgrp 
      use m_struc_def  !Cgetarg
      use m_hamindex,only: ibastab,ibasindex,symops,offl,offlrev,ltab,ktab,ag,invgx,nbas,
     &  shtvg,qq,qtt,qtti,iqmap,igmap,iqimap,dlmm,norbmto,ngrp,nqtt,nqi,plat,qlat,miat,tiat, !a little illegal
     &  ndimham,ngpmx,nqnum,iclasst,kxx,lxx,writehamindex,lxxa,
     &     napwk,napwmx,igv2,igv2rev, !for APW part of rotation (for rdsigm2.F)
     &     imx
     &  ,ngrpaf,symops_af,ag_af,iclasstaf, ngrp_original
      use m_lmfinit,only:
     &     ham_pwmode,ham_pwemax,
     &     lat_alat,nl,ctrl_nbas,ssite=>v_ssite,sspec=>v_sspec
      use m_lattic,only: lat_qlat,lat_plat,rv_a_opos
      use m_shortn3,only: shortn3_initialize,shortn3
!!-- Set up m_hamiltonian. Index for Hamiltonian. --
!!  Generated index are stored into m_hamindex
!!  Only include q-point information for GW (QGpsi).
!!
!!#Inputs
Co all output is in m_hamindex
C
Cr     iprmb is hard to understand (Mark's convension).  But, in anyway,
Cr     it is converted into clean indexing for Hamiltonian block.
Cr As you see in subroutine rotwv, the index for Hamiltonian reads as;
c      do iorb=1,norbmto             !orbital-blocks are specified by ibas, l, and k.
c        ibas  = ibastab(iorb)
c         l    = ltab(iorb)
c         k    = ktab(iorb)        !kappa (diffent MTO index for each l)
c        init1 = offl(iorb)+1      !starting index for the block iorb
c        iend1 = offl(iorb)+2*l+1  !end of the block for the iorb
c      enddo
!!h takao june2009
!!---
      implicit none
c      type(s_spec)::sspec(*)
c      type(s_site)::ssite(*)
c      type(s_lat)::slat
c     type(s_bz)::sbz
c      type(s_ham)::sham
c      integer iprmb(*),nl
      integer:: nkap0 , n0 , ibas , k 
      integer:: l , ndim , ipr , nglob ,  ldim , off , offs , 
     . specw , fieldw , iorb , offsi , ib , is , norb
      real(8):: alat,dum,qb(3,3),ddd(3) !plat(3,3),qlat(3,3),
      integer,allocatable:: ltabx(:,:),ktabx(:,:),offlx(:,:),iqtt(:)
      integer,parameter :: n0nkap0=30
      integer:: nkabc(3),nkp,lshft(3),pwmode,napwx,ig,nini
      integer:: nk1,nk2,nk3,ik1,ik2,ik3,ikt
      integer,allocatable:: kv(:)
      real(8):: pwgmax, pwemax, pwgmin, pwemin, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)
      integer:: ngp, ifiqg,iq,fopnx,nnn(3),ixx,ndummy,nqbz___ 
c      real(8):: qplist(3,nkp)
      logical:: siginit, qpgexist,debug=.false. !,no_qpgpsi
      integer:: i_copy_size,i_spacks,i_spackv,ifi,nbas_in,ifisym,
     &     i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw
      logical:: ljobgw

      character(8),allocatable::  spid(:)
      logical:: llmfgw
c      integer,allocatable:: igv2(:,:,:),napwk(:),igv2rev(:,:,:,:),ibasindex(:)
c      integer::  napwmx=null, nkt=null,ndimham=null
c      integer,allocatable:: ltab(:),ktab(:),offl(:), iclasst(:),offlrev(:,:,:),ibastab(:) !,ibasindex(:)
      integer::nkt
      logical:: iprx=.true.,prpushed

      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw ,iprint,ngadd,igadd,igaf
      real(8)::ppin(3)
      real(8):: tolq
      real(8):: rlatp(3,3),xmx2(3)
      real(8):: qqq(3),diffs,ddf
      real(8),allocatable:: symtmp(:,:,:)
      integer :: ifatomlist
#if (MPI|MPIK)
c2012Sep02 kino, add for 
       include 'mpif.h'
       integer:: procid=0,ier=0
       integer,parameter::master=0
       call mpi_comm_rank(mpi_comm_world,procid,ier)
       iprx=.false.
       if(procid==master) iprx=.true.
#endif
!!--- setup ---
      if(.not.siginit) return
      ldim  = ham_ldham(1)
      nbas=ctrl_nbas !nbas_in
      ngrp=lat_nsgrp !note nsgrp given in mksym.F is without inversion.
      plat=lat_plat
      qlat=lat_qlat
      alat=lat_alat
      allocate( ltabx(n0nkap0,nbas),ktabx(n0nkap0,nbas),offlx(n0nkap0,nbas)) !,ispec(nbas))
c      if(iprx) write(6,*) ' gen_hamindex: --- Hamiltonian index ---'
c      if(iprx) write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
!!--- MTO part ---
!!---- obtain norbmto, lxx,kxx ----
      norbmto=0
      kxx=-1
      lxx=-1
      ndimham = 0 !dimension of mto part of hamiltonian
      do  ib = 1, nbas
        is=ssite(ib)%spec
        call orbl(ib,0,ldim,iprmb,norb,ltabx(:,ib),ktabx(:,ib),off,offlx(:,ib),ndim)
        do iorb = 1, norb
          norbmto = norbmto+1
          if(ltabx(iorb,ib)>lxx)  lxx = ltabx(iorb,ib)
          if(ktabx(iorb,ib)>kxx)  kxx = ktabx(iorb,ib)
          ndimham = ndimham+ 2*ltabx(iorb,ib)+1
        enddo
      enddo
!!--- make index table ---
      allocate( ibasindex(ndimham))
      allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
      allocate(spid(nbas))
      norbmto=0
      ndimham = 0 !dimension of mto part of hamiltonian
      do  ib = 1, nbas
        is=ssite(ib)%spec
c        do i_spacks=is,is
c          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid(ib),i_spacks)
c        enddo
        spid(ib)=sspec(is)%name
        call orbl(ib,0,ldim,iprmb,norb,ltabx(:,ib),ktabx(:,ib),off,offlx(:,ib),ndim)
        do  iorb = 1, norb
          norbmto=norbmto+1
          ibastab(norbmto)= ib
          ltab(norbmto)   = ltabx(iorb,ib)
          ktab(norbmto)   = ktabx(iorb,ib)
          offl(norbmto)   = offlx(iorb,ib)
          nini = ndimham+ 1
          ndimham = ndimham+ 2*ltabx(iorb,ib)+1
          ibasindex(nini:ndimham) = ib
c          if(iprx)write(6,"(3i3,2x,2i5,3x,a)") 
c     &    ib,ltab(norbmto),ktab(norbmto), offl(norbmto)+1,offl(norbmto)+2*ltab(norbmto)+1,trim(spid)
        enddo
      enddo
      ifatomlist=999
      open(ifatomlist,file="ATOM_LIST")
      write(ifatomlist,*)nbas
      write(ifatomlist,*)spid
      close(ifatomlist)
c      print *, "SWJ", spid
c ... reverse maping of offset-index for hamiltonian
      allocate(offlrev(nbas,0:lxx,kxx))
      do iorb=1,norbmto
        ibas = ibastab(iorb)
        l   = ltab(iorb)
        k   = ktab(iorb)
        offlrev(ibas,l,k)= offl(iorb)
      enddo


!!--- symmetry operation ---
      allocate(symops(3,3,ngrp),ag(3,ngrp))
      call dcopy ( ngrp * 9 , rv_a_osymgr , 1 , symops , 1 )
      call dcopy ( ngrp * 3 , rv_a_oag , 1 , ag , 1 )
      allocate(iclasst(nbas),invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      do ib=1,nbas
        iclasst(ib)=ssite(ib)%class
      enddo
!!---- get space group information ---- translation informations also in miat tiat invgx, shtvg
      call mptauof ( symops , ngrp , plat , nbas , rv_a_opos , iclasst 
     . , miat , tiat , invgx , shtvg )
!!---- write SYMOPS ---- mar2012takao
#if (MPI|MPIK)
      if (procid.eq.master) then
#endif
c     ifisym = fopnx('SYMOPS',2,0,-1) ! open as unknown
      open(newunit=ifisym,file='SYMOPS')
      write(ifisym,*) ngrp
      do ig = 1,ngrp
        write(ifisym,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
        do i=1,3
          write(ifisym,"(3e24.16)") symops(i,1:3,ig)
        enddo
      enddo
      close(ifisym)
#if (MPI|MPIK)
      endif
#endif

!! Get rotation matrix dlmm.  We assume nl=lmxa+1.
      lxxa=nl-1
      allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      call rotdlmm(symops,ngrp, nl, dlmm)

!! No GW-related part
      if(ljobgw) then
c      if((.not.llmfgw).or.(llmfgw.and.jobgw==0)) then
        if(iprx) print *,'gen_hamindex: not readin QGpsi.'
        return
      endif
!! --- PW part. info for eigenfunctions are expanded as MTpart+PWpart.!feb2012takao
      inquire(file='QGpsi',EXIST=qpgexist)  !feb2012takao
      if(.not.qpgexist) then
        call writehamindex()
        goto 2001 !return !jun2015takao
      endif
!! q on mesh and shortened q.
      open(newunit=ifiqg,file='QGpsi',form='unformatted')
      read(ifiqg) nqnum, ngpmx ,QpGcut_psi, nqbz___, nqi !,imx !,nqibz
c$$$cccccccccccccccccccccccccccccc
c$$$      nkt = 2*nqbz + 2*nkp
c$$$      if(allocated(qq)) deallocate(qq)
c$$$      allocate( qq(3,nkt) )
c$$$      print *,'gen_hamindex: nkt nqbz=',nkt,nqbz
c$$$      do  iq = 1, nqbz
c$$$        read(ifiqg)  qxx  ! q, and number of G vectors for
c$$$        read(ifiqg)
c$$$        qq(:,iq)=qxx
c$$$        call shorbz(qq(:,iq),qq(:,iq+nqbz+nkp),qlat,plat)
c$$$        write(6,"(' qq=',i5,3f10.5,' shortened qq=',3f10.5)") iq,qq(:,iq),qq(:,iq+nqbz+nkp)
c$$$      enddo
c$$$      print *
c$$$      call fclr(' ',ifiqg)
c$$$      do iq = nqbz+1, nqbz+nkp
c$$$        qq(:,iq) = qplist(:,iq-nqbz)
c$$$        call shorbz(qq(:,iq),qq(:,iq+nqbz+nkp),qlat,plat)
c$$$        write(6,"(' qq=',i5,3f10.5,' shortened qq=',3f10.5)") iq,qq(:,iq),qq(:,iq+nqbz+nkp)
c$$$      enddo
cccccccccccccccccccccccccccccc

!! === feb2012takao. ===
!! we have two set of data for original qxx in QGpsi and their shortened.
c      nqnum2 = 2*nqnum !+ 2*nkp !takao feb2012 test xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(allocated(qq)) deallocate(qq)
!! aug2012 we are removing shorbz
      nqtt=nqnum !nqnum*2 !doubled. second series nqnum+1:2*nqnum are for shortened q.
      nkt=nqtt
      allocate( qtti(3,nqi), qq(3,nqtt),iqtt(nqtt) )
c      allocate( ngvecp(3,ngpmx,nqnum))
c      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
      if(iprx) print *
      if(iprx) print *,'gen_hamindex: Readin QGpsi. nqnum=',nqnum
      iqi=0
      do  iq = 1, nqnum
        read(ifiqg)  qxx,ngp,irr  ! q, and number of G vectors for
        if(irr/=0) then
          iqi=iqi+1
          qtti(:,iqi)=qxx
          iqtt(iqi)=iq
        endif
        read(ifiqg)
c        read(ifiqg) ngvecp(1:3,1:ngp,iq), ngvecprev(-imx:imx,-imx:imx,-imx:imx,iq)
        qq(:,iq)=qxx
cc comment out aug2012takao
cc        call shorbz(qq(:,iq),qq(:,iq+nqnum),qlat,plat)
cc        write(6,"(' qq=',i5,3f10.5,' shortened qq=',3f10.5)") iq,qq(:,iq),qq(:,iq+nqnum)
        if(iprx)write(6,"(' qq=',i5,3f10.5)") iq,qq(:,iq)
      enddo
      close(ifiqg)
cccccccccccccccccc
c        do iqi=1,nqi
c          q=qtti(:,iqi)
c          print *,' qtti=',iqi,q
c        enddo
c        stop 'xxxxxxxxxxxxxxxx'
cccccccccccccccccc

!! ==== Generate info for rotwv and write ====
c     print *,' nqtt nqi=',nqtt,nqi
c     print *, ' symops=',symops
      allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
      platt= transpose(plat) !this is inverse of qlat
      allocate(qtt(3,nqtt))
      qtt(:,1:nqtt)=qq(:,1:nqtt)
      do i=1,nqtt
        qtarget(:)=qtt(:,i)
c        print *,' i in nqtt=',i,qtarget
        do iqi=1,nqi
          q=qtti(:,iqi)
          iqq=iqtt(iqi)
          iqi_=iqi
c          print *,' xxxxx qtti=',q
          do ig=1,ngrp
c       print *, 'ddd=',matmul(platt,(qtarget-matmul(symops(:,:,ig),q)))
c            qx=matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) )
c            if(sum(abs(qx-anint(qx)))<1d-8) then
            call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
            if(sum(abs(qx))<tolq()) then
              igg=ig
c              if(debug)print *,' q qtarget     =',q,qtarget,ig
c              if(debug)print *,' matmul q =',matmul(symops(:,:,ig),q)
c              if(debug)print *
              ddf=sum(abs(matmul(platt,(qtarget-matmul(symops(:,:,ig),q)))))
              if(ddf-nint(ddf) > 1d-8) then
              write(6,"('qqqqqxq ',3d16.8,2x,3d16.8,2x,3d16.8)") q, qtarget,qtarget-matmul(symops(:,:,ig),q)
              endif
              goto 2012
            endif
          enddo
        enddo

        if(iprx) then
           write(6,"(a,3f7.3,2x,3f7.3)")'gen_ham: qtarget cannot found.'//
     &     ' Need to add SYMGRP explicitly (for SO=1), or You have to delete inconsistent QGpsi. qtarget=',qtarget
           print *,'gen_hamindex: qtarget can not found by SYMOPS.' 
           write(6,"('qqqqqxq20 ',3d16.8,2x,3d16.8,2x,3d16.8)") q, qtarget, matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) )
           do ig=1,ngrp
           call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
           write(6,"('qqqqqxq2 ',3d16.8,2x,3d16.8,2x,3d16.8)") qtarget-matmul(symops(:,:,ig),q),qx
           enddo
           call rx('gen_hamindex: you may need to repeat echo 1|qg4gw, when you changed SYMOPS.')
        endif

 2012   continue
        iqmap(i)=iqq
        iqimap(i)=iqi_
        igmap(i)=igg
      enddo

!! === rotation of APW. (not the IPW part for GW).===
c$$$ 1012 continue
      pwmode=ham_pwmode
      pwemax=ham_pwemax
c     pwemin=ham_pwemin
      if(pwmode==0.or.pwemax<1d-8) then
        if(allocated(napwk)) deallocate(napwk)
        allocate(napwk(nkt))
        napwk=0
        napwmx=0
        print *,'pwmode=0 writehamindex'
        call writehamindex() !sep2012takao
        return
      endif

!! for APW rotation.
C ... Get igv2(3,iapw,ikt). pwmode>=10 only
      if(iprx) print *,' gen_hamindex goto APW part: pwmode pwemax=',pwmode,pwemax !pwemin
      if(allocated(napwk)) deallocate(napwk,igv2,igv2rev)
      allocate( napwk(nkt))
!! takao is 
      if(mod(pwmode,10)==0) then ! MTO basis only 
        return
      endif
      pwgmax = dsqrt(pwemax)
      pwgmin = 0d0 !dsqrt(pwemin) !this will be removed.
      napwmx = 0
      call pshpr(0) !print index is pushed to be zero 
      do ikt=1,nkt
        qqq=0d0 !call dpzero(xx,3)
        if (mod(pwmode/10,10) == 1) qqq=qq(:,ikt) !call dpcopy(qp,xx,1,3,1d0)
        call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,0,napwx,dum,dum,dum,dum)
        napwk(ikt) = napwx
        if(napwmx<napwx) napwmx = napwx
      enddo
      call poppr
cnnn
      if(iprx)print*,' --- gvlst2 generates G for APW part (we show cases for limited q) ---'
      if(pwmode<5) call shortn3_initialize(qlat)
      allocate( igv2(3,napwmx,nkt), kv(3*napwmx) )
      prpushed=.false.
      do ikt = 1,nkt
        qqq=0d0
        if (mod(pwmode/10,10) == 1) qqq=qq(:,ikt) 
        call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napwmx,napwk(ikt),kv,dum,dum,igv2(1,1,ikt))
        if(iprx.and.(ikt>5.or.ikt==nkt).and.(.not.prpushed)) then
          call pshpr(0)
          prpushed=.true.
        endif
cnnn
        if(pwmode<10) then
        ppin=matmul(transpose(plat),qq(:,ikt))
        call shortn3(ppin,noutmx, nout,nlatout)
c        do iout=1,nout
c           write(*,"(a,3i5,f10.4,3f8.4)")'suham.gen_hamindex.F:rrrrn1 =',nlatout(:,iout),
c     &         sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
c     &         matmul(qlat(:,:),ppin+nlatout(:,iout))
c        enddo
        endif
cnnn
        if (pwmode<10) then
           do iapw=1,napwk(ikt)
              igv2(:,iapw,ikt)=igv2(:,iapw,ikt)+nlatout(:,1)
           enddo
           
        endif
      enddo
      deallocate(kv)
      if(iprx) call poppr !print index is poped.

C ... Reverse table of igv2 --->igv2rev
      imx=-999
      do ikt = 1,nkt
        ixx = maxval( abs(igv2(1:3,1:napwk(ikt),ikt)))
        if(ixx>imx) imx=ixx
      enddo
      allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt) )
      igv2rev=999999
      do ikt = 1,nkt
cccccccccccccccccccccccc
c          write(6,"(a,3f8.3,i4)")'rrrrrrr suham.gen_hamindex qq ikt: ',qq(:,ikt),ikt
cccccccccccccccccccccccc
        do ig  = 1,napwk( ikt )
          nnn  = igv2(1:3, ig, ikt)
          igv2rev( nnn(1), nnn(2),nnn(3), ikt) = ig
ccccccccccccccccccccccc
c          write(6,"(a,3f8.3,4i4,i6,i6)")'suham.gen_hamindex: igv2rev: ',qq(:,ikt), nnn, ig, ikt
c          write(6,"(a,i4,2x,3i4)")'suham.gen_hamindex: igv2: ', ig, nnn
ccccccccccccccccccccccc
        enddo
      enddo
      call writehamindex() !sep2012takao
      if(iprx) print *,' end of gen_hamindex: nkt,napwmx norbmto= ',nkt,napwmx,norbmto

!! Symmetry for AF. Order AF symmetry operation after normal one. jun2015takao
!! Caution: symops,and so on are overwritten.
!!          writehamindex() already wrote HAMindex which is just for SYMGRP.
 2001 continue
      if(allocated(symops_af)) then
        allocate(symtmp(3,3,ngrpaf))
        symtmp(:,:,1:ngrp)=symops
        igadd=ngrp
        do igaf=1,ngrpaf
        do ig=1,ngrp
          diffs=sum(abs(symops_af(:,:,igaf)-symops(:,:,ig)))
          if(diffs<1d-6) then
            goto 1013
          endif  
        enddo
        igadd=igadd+1
c       if(iprx)  write(6,*) ' add igaf=',igaf
        symtmp(:,:,igadd)=symops_af(:,:,igaf)
 1013   continue
        enddo
        if(igadd/=ngrpaf) call rx('suham: strange. bug igadd/=ngrpaf')
        if(iprx) write(6,*) '-----SYMGRPAF mode ---- # of additional symmetry=',igadd
!!---- get space group information ---- translation informations also in miat tiat invgx, shtvg
        deallocate(symops, invgx,miat,tiat,shtvg,  dlmm )
        ngrp_original=ngrp
        ngrp      = ngrpaf ! Overwrite ngrp nby ngrpaf (larger, because we assume AF site are in
                           ! the same class).
        allocate(symops(3,3,ngrp))
        symops = symtmp
        allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
        call mptauof ( symops , ngrp, plat , nbas , rv_a_opos , iclasstaf
     .   , miat , tiat , invgx , shtvg )
        if(iprx) then
          write(6,*)
          write(6,"(' ngrp for SYMGRP+GYMGRPAF, ngrp for SYMGRP=',2i5)") ngrp, ngrp_original
          do ig=1,ngrp
          write(6,"(a,i3,a,100i3)")'ig=',ig, ' miat=',miat(:,ig)
          enddo
        endif  
!! Get rotation matrix dlmm.  We assume nl=lmxa+1.
        lxxa=nl-1
        allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
        call rotdlmm(symops,ngrp, nl, dlmm)
c        call rx0('xxxxxxxxxxxx test end af mode for suham xxxxxxxxxxxxx')
      endif
      end subroutine gen_hamindex
csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine mptauof(symops,ng,plat,nbas,bas,
     &iclass,miat,tiat,invg,delta)

C- Mapping of each atom by points group operations------------------c
C Modded by okuda 1994 March.
C Simplified by kotani 1994 7/31.
Ci  Input
Ci     symops(1,ng),ng,plat,nbas,bas(3,nbas)
Ci     iclass(nbas); denote class for each atom
Co  Output
Co    miat(ibas  ,ig); ibas-th atom is mapped to miat-th atom, by the ig-th
Co    points group operation.  Origin is (0,0,0).
Co    tiat(k,ibas,ig);
Co    delta : shifting vector for non-symmorphic group.
Co            r' = matmul (am, r) + delta
Cr  Remarks
Cr
Cr (1) The ibas-th atom (position at bas(k,ibas) ) is mapped to
Cr
Cr    bas( k,miat(ibas,ig) )+ tiat(k,ibas,ig), k=1~3.
Cr
Cr (2) tiat= unit translation
Cr
C--------------------------------------------------------------------
      implicit none
      integer ng,nbas, miat(nbas,ng),iclass(nbas),invg(ng),
     &nbmx, nsymx, ig,igd,i,j,ibas,mi,i1,i2,i3
      double precision SYMOPS(9,ng),plat(3,3),
     &tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas),
     &tr1,tr2,tr3,ep, dd1,dd2,dd3,t1,t2,t3
      integer::  iprintx=0
      integer ires(3, nbas, ng)
ckino delete integer(4) def.      integer(4):: ib1,ib2
      integer:: ib1,ib2
c
      real(8) ::tran(3),delta(3,ng)
      data ep/1.0d-3/
c      data ep/1.0d-7/
c
      if(iprintx>=46) write(6,*)'MPTAUOf: search miat tiat for wave function rotation'

      do 10 ig=1,ng
        do igd=1,ng
c seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ).le.ep.and.
     &    abs( symops(2,ig)-symops(4,igd) ).le.ep.and.
     &    abs( symops(3,ig)-symops(7,igd) ).le.ep.and.
     &    abs( symops(4,ig)-symops(2,igd) ).le.ep.and.
     &    abs( symops(5,ig)-symops(5,igd) ).le.ep.and.
     &    abs( symops(6,ig)-symops(8,igd) ).le.ep.and.
     &    abs( symops(7,ig)-symops(3,igd) ).le.ep.and.
     &    abs( symops(8,ig)-symops(6,igd) ).le.ep.and.
     &    abs( symops(9,ig)-symops(9,igd) ).le.ep  ) then
            invg(ig)=igd
            goto 16
          endif
        end do
 16     continue
c
        if(iprintx .ge.46) then
          print *,' '
          print *,' '
          print *,' **** group ops no. ig (igd)= ', ig, invg(ig)
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
 1731     format (' ',3f9.4)
        endif

        do i=1,3
          do j=1,3
            am(i,j)=symops(i+3*(j-1),ig)
          end do
        end do
c
c trial shift vector tran
        do 20 ib1=1,nbas
        do 20 ib2=1,nbas
          tran =  bas(:,ib2)  - matmul(am,bas(:,ib1))
c
          do 30 ibas=1,nbas
            b1=am(1,1)*bas(1,ibas)
     &        +am(1,2)*bas(2,ibas)+am(1,3)*bas(3,ibas)
     &        +tran(1)
c     .        +( tr1*plat(1,1)+tr2*plat(1,2)+tr3*plat(1,3) )
            b2=am(2,1)*bas(1,ibas)
     &        +am(2,2)*bas(2,ibas)+am(2,3)*bas(3,ibas)
     &        +tran(2)
c     .        +( tr1*plat(2,1)+tr2*plat(2,2)+tr3*plat(2,3) )
            b3=am(3,1)*bas(1,ibas)
     &        +am(3,2)*bas(2,ibas)+am(3,3)*bas(3,ibas)
     &        +tran(3)
c     .        +( tr1*plat(3,1)+tr2*plat(3,2)+tr3*plat(3,3) )
c
            do 40 mi=1,nbas
              if( iclass(mi).ne.iclass(ibas) ) go to 40

              do 50 i1=-3,3
              do 50 i2=-3,3
              do 50 i3=-3,3
                dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )

                t1 = b1 - (bas(1,mi)+dd1)
                t2 = b2 - (bas(2,mi)+dd2)
                t3 = b3 - (bas(3,mi)+dd3)
                if(abs(t1).le.ep.and.abs(t2).le.ep.and.
     &                abs(t3).le.ep) go to 60
   50         continue
   40       continue
c seach failed, Not found mi and dd1. Try next (tr).
            goto 20

   60       continue
            miat(ibas,ig)  = mi
            tiat(1,ibas,ig)= dd1
            tiat(2,ibas,ig)= dd2
            tiat(3,ibas,ig)= dd3
            ires(1,ibas,ig)= i1
            ires(2,ibas,ig)= i2
            ires(3,ibas,ig)= i3
c
   30     continue
c When the do-30 loop has been completed, we get out of do-20 loop
          goto 21
   20   continue
        call rx('mptauof: Can not find miat and tiat')
c
   21   continue
        delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000

c- have gotten the translation-> check write --------------------
        if(iprintx.ge.46) then
c          write(6,4658)tr1,tr2,tr3
          write(6,4658)tran
 4658     format('  Obtained translation operation=',3d12.4)
          do 123  ibas=1,nbas
            write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig),
     &      tiat(2,ibas,ig), tiat(3,ibas,ig),
     &      ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
  150       format(' ibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
  123     continue
        endif
c---------------------------------------------------
   10 continue
      end subroutine mptauof


      subroutine rotdlmm(symops,ng,nl ,dlmm)
c- Generate rotation matrix D^l_{m,m'} for L-representaiton, corresponding
c  to points group operations.
Ci symops(9,ng),ng; point ops.
Ci nl; num.of l =lmax+1
Co dlmm(2*nl-1,2*nl-1,0:nl-1,ng,2); D^l_{m,m'}. Indexes are for Real harmonics.
Cr   dlmmc is used as work area about 200kbyte used for  s,p,d,f -> nl=4
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer:: is,i,ig,ikap,j,l,m,m1,m2,m3,md,mx,ix
      integer ng,nl
      double precision SYMOPS(9,ng),
     .am(3,3) ,fac1,fac2
      double precision dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)

      double complex   dlmmc(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
      double precision det,igann,osq2
      double complex   msc(0:1,2,2), mcs(0:1,2,2),Img
     &,dum(2)
      parameter(Img=(0d0,1d0))
      integer:: debugmode
      real(8):: ep=1d-3 !ep was 1d-8 before feb2013
c      print *; print *,' ROTDLMM:'
      do 10 ig =1,ng
        do 20 i=1,3
        do 20 j=1,3
          am(i,j) = symops(i+3*(j-1),ig)
   20   continue
c calculate determinant(signature)
        det= am(1,1)*am(2,2)*am(3,3)
     &  -am(1,1)*am(3,2)*am(2,3)
     &  -am(2,1)*am(1,2)*am(3,3)
     &  +am(2,1)*am(3,2)*am(1,3)
     &  +am(3,1)*am(1,2)*am(2,3)
     &  -am(3,1)*am(2,2)*am(1,3)
        if(abs(abs(det)-1d0).ge.1d-10) then
          print *,' rotdlmm: det/=1 ig and det=',ig,det
          stop
        endif
c seek Euler angle   print *,' goto cbeta',ig,det
        cbeta = am(3,3)/det
c added region correction so as to go beyond domain error for functions, dsqrt and acos.
        if(abs(cbeta-1d0).le.1d-6) cbeta= 1d0
        if(abs(cbeta+1d0).le.1d-6) cbeta=-1d0
        beta = dacos(cbeta)
        sbeta= sin(beta)
c beta= 0~pi
        if(sbeta.le.1.0d-6) then
          calpha= 1d0
          salpha= 0d0
          alpha = 0d0
          cgamma= am(2,2)/det
          sgamma= am(2,1)/det
        else
          salpha =  am(2,3)/sbeta/det
          calpha =  am(1,3)/sbeta/det
          sgamma =  am(3,2)/sbeta/det
          cgamma = -am(3,1)/sbeta/det
        endif
        co2 = dcos(beta/2d0)
        so2 = dsin(beta/2d0)
c         print *,' calpha=',calpha
        if(abs(calpha-1.0d0).le.1.0d-6) calpha= 1.0d0
        if(abs(calpha+1.0d0).le.1.0d-6) calpha=-1.0d0
        if(abs(cgamma-1.0d0).le.1.0d-6) cgamma= 1.0d0
        if(abs(cgamma+1.0d0).le.1.0d-6) cgamma=-1.0d0
        alpha=dacos(calpha)
        if(salpha.lt.0d0) alpha=-alpha
        gamma=dacos(cgamma)
        if(sgamma.lt.0d0) gamma=-gamma
c         print *,'alpha beta gamma det=',alpha,beta,gamma,det
        do 30 l =  0, nl-1
        do 30 md= -l, l
        do 30 m = -l, l
c  from 'Ele theo. ang. mom. by M. E. Rose 5th 1967 Wisley and Sons.  p.52 (4.13)
          fac1 = dsqrt( igann(l+m)*igann(l-m)*igann(l+md)*igann(l-md) )
          fac2 = 0d0
          do 40 ikap=0,2*l
            if(l-md-ikap.ge.0 .and. l+m-ikap.ge.0
     &          .and.ikap+md-m.ge.0) then
              add= dble((-1)**ikap)/( igann(l-md-ikap)*igann(l+m-ikap)
     &            *igann(ikap+md-m)*igann(ikap) )
              if(2*l+m-md-2*ikap.ne.0) add=add*co2**(2*l+m-md-2*ikap)
              if(md-m+2*ikap.ne.0)     add=add*(-so2)**(md-m+2*ikap)
              fac2 = fac2+add
            endif
   40     continue
c l-th rep. is odd or even according to (det)**l
          dlmmc(md,m,l,ig) = fac1*fac2*det**l*
     &        cdexp( -Img*(alpha*md+gamma*m) )
   30   continue

        am(1,1)= cos(beta)*cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)
        am(1,2)=-cos(beta)*cos(alpha)*sin(gamma)-sin(alpha)*cos(gamma)
        am(1,3)= sin(beta)*cos(alpha)
        am(2,1)= cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
        am(2,2)=-cos(beta)*sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)
        am(2,3)= sin(beta)*sin(alpha)
        am(3,1)=-sin(beta)*cos(gamma)
        am(3,2)= sin(beta)*sin(gamma)
        am(3,3)= cos(beta)

        if(abs(am(1,1)*det-symops(1,ig))>ep.or.
     &  abs(am(2,1)*det-symops(2,ig))>ep.or.
     &  abs(am(3,1)*det-symops(3,ig))>ep.or.
     &  abs(am(1,2)*det-symops(4,ig))>ep.or.
     &  abs(am(2,2)*det-symops(5,ig))>ep.or.
     &  abs(am(3,2)*det-symops(6,ig))>ep.or.
     &  abs(am(1,3)*det-symops(7,ig))>ep.or.
     &  abs(am(2,3)*det-symops(8,ig))>ep.or.
     &  abs(am(3,3)*det-symops(9,ig))>ep) then
          print *,' rotdlmm: not agree. symgrp and one by eular angle'
          stop
        endif
cccccccccccccccccccccccc
c        if(iprint().ge.140) then
        if(debugmode()>9) then
          print *;print *;print *,' **** group ops no. ig=', ig
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
          print *,' by Eular angle '
          write(6,1731)am(1,1)*det,am(1,2)*det,am(1,3)*det
          write(6,1731)am(2,1)*det,am(2,2)*det,am(2,3)*det
          write(6,1731)am(3,1)*det,am(3,2)*det,am(3,3)*det
        endif
 1731   format (' ',3f9.4)
cccccccccccccccccccccccc
   10 continue
c conversion to cubic rep. Belows are from csconvs
c  msc mcs conversion matrix generation 2->m 1->-m for m>0
      osq2 = 1d0/sqrt(2d0)
      do m = 0,1
        Msc(m,1,1)= osq2*(-1)**m
        Msc(m,1,2)=-osq2*Img*(-1)**m
        Msc(m,2,1)= osq2
        Msc(m,2,2)= osq2*Img

        Mcs(m,1,1)= osq2*(-1)**m
        Mcs(m,1,2)= osq2
        Mcs(m,2,1)= osq2*Img*(-1)**m
        Mcs(m,2,2)=-osq2*Img
      enddo
c
      if(debugmode()>1) print * ,' goto do 23'
      do 23 is=1,ng
        if(.false.) then
c        if(iprint().ge.150) then
          print *; print *,' **** group ops no. ig=', is
          write(6,1731) symops(1,is),symops(4,is),symops(7,is)
          write(6,1731) symops(2,is),symops(5,is),symops(8,is)
          write(6,1731) symops(3,is),symops(6,is),symops(9,is)
        endif
c convert to cubic rep.
      do 23   l =0,nl-1
        do 33 m2=-l,l
        do 33 m1= 1,l
          dum(1)= dlmmc(m2, m1,l,is)
          dum(2)= dlmmc(m2,-m1,l,is)
          mx    = mod(m1,2)
          dlmmc(m2,  m1,l,is)=
     &        dum(1)*msc(mx,1,1)
     &        +dum(2)*msc(mx,2,1)
          dlmmc(m2, -m1,l,is)=
     &        dum(1)*msc(mx,1,2)
     &        +dum(2)*msc(mx,2,2)
   33   continue
        do 43 m2=  1,l
        do 43 m1= -l,l
          dum(1)=dlmmc( m2, m1,l,is)
          dum(2)=dlmmc(-m2, m1,l,is)
          mx=mod(m2,2)
          dlmmc( m2, m1,l,is)=
     &        mcs(mx,1,1)*dum(1)
     &        +mcs(mx,1,2)*dum(2)
          dlmmc(-m2, m1,l,is)=
     &        mcs(mx,2,1)*dum(1)
     &        +mcs(mx,2,2)*dum(2)
   43   continue
        do 53 m2=-l,l
        do 53 m1=-l,l
          dlmm(m2,m1,l,is)=dreal( dlmmc(m2,m1,l,is) )
          if( abs(dimag(dlmmc(m2,m1,l,is))).ge.1.0d-12 ) stop
     &        ' rotdlmm: abs(dimag(dlmmc(m2,m1,l,is))).ge.1.0d-12'
   53   continue
ccccccccccccccccccccc
        if(.false.) then
c        if(.true.) then
c        if(iprint().ge.41) then
          print *; print *,'  points ops  ig, l=', is,l,' cubic   '
          do m2=-l,l
            write(6,"(28f10.5)")( dreal(dlmmc (m2, m1,l,is) ), m1=-l,l)
c    &    , ( dimag(dlmmc (m2, m1,l,is) ), m1=-l,l),( dlmm(m2, m1,l,is), m1=-l,l)
          enddo
        endif
cccccccccccccccccccccc
   23 continue
      if(debugmode()>1) print *,' end of rotdlmm'
      end subroutine rotdlmm

c--------------------------------------------
      double precision function igann(i)
      integer:: i,ix
      igann  = 1d0
      do ix =1,i
        igann=igann*dble(ix)
      enddo
      end function igann

