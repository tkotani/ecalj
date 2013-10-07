!! we have to keep these( m_hamiltoninan ) with those in fpgw/gwsrc
      subroutine gen_hamindex(sham,sbz,ssite,sspec,slat,nl,nbas_in,iprmb,ldim,siginit,llmfgw,jobgw) !,no_qpgpsi) !from showbs
c      subroutine gen_hamindex(sham,sbz,ssite,sspec,slat,nl,nbas_in,iprmb,ldim,qplist,nkp,siginit) !,no_qpgpsi) !from showbs
      use m_struc_def  !Cgetarg
      use m_hamindex,only: ibastab,ibasindex,symops,offl,offlrev,ltab,ktab,ag,invgx,nbas,
     &  shtvg,qq,qtt,qtti,iqmap,igmap,iqimap,dlmm,norbmto,ngrp,nqtt,nqi,plat,qlat,miat,tiat,
     &  imx,ndimham,ngpmx,nqnum,iclasst,kxx,lxx,writehamindex,lxxa,
     &  napwk,napwmx,igv2,igv2rev !for APW part of rotation (for rdsigm2.F)

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
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)
      type(s_lat)::slat
      type(s_bz)::sbz
      type(s_ham)::sham
      integer iprmb(*),nl
      integer:: nkap0 , n0 , ibas , k
      integer:: l , ndim , ipr , nglob ,  ldim , off , offs , 
     . specw , fieldw , iorb , offsi , ib , is , norb
      character spid*8
      real(8):: alat,dum,qb(3,3),ddd(3) !plat(3,3),qlat(3,3),
      integer,allocatable:: ltabx(:,:),ktabx(:,:),offlx(:,:),iqtt(:)
      integer,parameter :: n0nkap0=30
      integer:: nkabc(3),nkp,lshft(3),pwmode,napwx,ig,nini
      integer:: nk1,nk2,nk3,ik1,ik2,ik3,ikt
      integer,allocatable:: kv(:)
      real(8):: pwgmax, pwemax, pwgmin, pwemin, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)
      integer:: ngp, ifiqg,iq,fopnx,nnn(3),ixx,ndummy,nqbz___ !imx,
c      real(8):: qplist(3,nkp)
      logical:: siginit, qpgexist,debug=.false. !,no_qpgpsi
      integer:: i_copy_size,i_spacks,i_spackv,ifi,nbas_in,ifisym,
     &  i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw

      logical:: llmfgw
c      integer,allocatable:: igv2(:,:,:),napwk(:),igv2rev(:,:,:,:),ibasindex(:)
c      integer::  napwmx=null, nkt=null,ndimham=null
c      integer,allocatable:: ltab(:),ktab(:),offl(:), iclasst(:),offlrev(:,:,:),ibastab(:) !,ibasindex(:)
      integer::nkt
      logical:: iprx=.true.

      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw
      real(8)::ppin(3)
      real(8):: rlatp(3,3),xmx2(3)
      real(8):: qqq(3)
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
      nbas=nbas_in
      ngrp=slat%nsgrp
      plat=slat%plat
      qlat=slat%qlat
      alat=slat%alat
      allocate( ltabx(n0nkap0,nbas),ktabx(n0nkap0,nbas),offlx(n0nkap0,nbas)) !,ispec(nbas))
      if(iprx) write(6,*) ' gen_hamindex: --- Hamiltonian index ---'
      if(iprx) write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
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
      norbmto=0
      ndimham = 0 !dimension of mto part of hamiltonian
      do  ib = 1, nbas
        is=ssite(ib)%spec
        do i_spacks=is,is
          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
        enddo
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
          if(iprx)write(6,"(3i3,2x,2i5,3x,a)") 
     &    ib,ltab(norbmto),ktab(norbmto), offl(norbmto)+1,offl(norbmto)+2*ltab(norbmto)+1,trim(spid)
        enddo
      enddo
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
      call dcopy ( ngrp * 9 , slat%rv_a_osymgr , 1 , symops , 1 )
      call dcopy ( ngrp * 3 , slat%rv_a_oag , 1 , ag , 1 )
      allocate(iclasst(nbas),invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      do ib=1,nbas
        iclasst(ib)=ssite(ib)%class
      enddo
!!---- get space group information ----
!!# translation informations also in miat tiat invgx, shtvg
      call mptauof ( symops , ngrp , plat , nbas , slat%rv_a_opos , iclasst 
     . , miat , tiat , invgx , shtvg )

!!---- write SYMOPS ---- mar2012takao
#if (MPI|MPIK)
      if (procid.eq.master) then
#endif
ckino2012Sep02      ifisym = fopnx('SYMOPS',2,2,-1) ! open as new, it causes a stop.
      ifisym = fopnx('SYMOPS',2,0,-1) ! open as unknown
      write(ifisym,*) ngrp
      do ig = 1,ngrp
        write(ifisym,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
        do i=1,3
          write(ifisym,"(3e24.16)") symops(i,1:3,ig)
        enddo
      enddo
ckino2012Sep02      call fclr(' ',ifisym)
      call fclose(ifisym)
#if (MPI|MPIK)
      endif
#endif
c ... get rotation matrix dlmm.  We assume nl=lmxa+1.
      lxxa=nl-1
      allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      call rotdlmm(symops,ngrp, nl, dlmm)
c      print *,' end of rotdlmm'
!! no GW-related part
c      if(.not.llmfgw) then !.or.(llmfgw.and.jobgw==0)) then
      if((.not.llmfgw).or.(llmfgw.and.jobgw==0)) then
        if(iprx) print *,'gen_hamindex: not readin QGpsi.'
        return
      endif
!! --- PW part. info for eigenfunctions are expanded as MTpart+PWpart.!feb2012takao
      inquire(file='QGpsi',EXIST=qpgexist)  !feb2012takao
      if(.not.qpgexist) then
        return
      endif
!! q on mesh and shortened q.
      ifiqg  = fopnx('QGpsi',2,4,-1)
      read(ifiqg) nqnum, ngpmx ,QpGcut_psi, nqbz___, nqi,imx !,nqibz
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
      call fclr(' ',ifiqg) !close file
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
            call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
            if(sum(abs(qx))<1d-8) then
              igg=ig
              if(debug)print *,' q qtarget     =',q,qtarget,ig
              if(debug)print *,' matmul q =',matmul(symops(:,:,ig),q)
              if(debug)print *
              goto 2012
            endif
          enddo
        enddo
        if(iprx)write(6,"(a,3f7.3,2x,3f7.3)")'gen_ham: qtarget cannot found.'//
     & ' Need to add SYMGRP explicitly (for SO=1), or You have to delete inconsistent QGpsi. qtarget=',qtarget

        print *,'gen_hamindex: qtarget can not found by SYMOPS.' 
        call rx('gen_hamindex: you may need to repeat echo 1|qg4gw, when you changed SYMOPS.')
 2012   continue
        iqmap(i)=iqq
        iqimap(i)=iqi_
        igmap(i)=igg
      enddo

!! === rotation of APW. (not the IPW part for GW).===
c$$$ 1012 continue
      pwmode=sham%pwmode
      pwemax=sham%pwemax
c     pwemin=sham%pwemin
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
      call pshpr(0)
c      print *,' zzzzzzzzzzz nkt=',nkt
      do ikt=1,nkt
c        print *,' 111 ikt nkt=',ikt,nkt,qq(:,ikt)
        !print *,' 111 alat plat=',alat,plat
        !print *,' 111 qq pwgmin pwgmax =',qq(:,ikt),pwgmin,pwgmax
        qqq=0d0 !call dpzero(xx,3)
        if (mod(pwmode/10,10) == 1) qqq=qq(:,ikt) !call dpcopy(qp,xx,1,3,1d0)
        call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,0,napwx,dum,dum,dum,dum)
        napwk(ikt) = napwx
        if(napwmx<napwx) napwmx = napwx
      enddo
      call poppr
cnnn
      if(pwmode<10) call shortn3_initialize(qlat,rlatp,xmx2)
      allocate( igv2(3,napwmx,nkt), kv(3*napwmx) )
      do ikt = 1,nkt
        qqq=0d0 !call dpzero(xx,3)
        if (mod(pwmode/10,10) == 1) qqq=qq(:,ikt) !call dpcopy(qp,xx,1,3,1d0)
        call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napwmx,napwk(ikt),kv,dum,dum,igv2(1,1,ikt))
cnnn
        if(pwmode<10) then
        ppin=matmul(transpose(plat),qq(:,ikt))
        call shortn3(rlatp,xmx2,ppin,noutmx, nout,nlatout)
        do iout=1,nout
           write(*,"(a,3i5,f10.4,3f8.4)")'suham.gen_hamindex.F:rrrrn1 =',nlatout(:,iout),
     &         sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
     &         matmul(qlat(:,:),ppin+nlatout(:,iout))
        enddo
        endif
cnnn
        if (pwmode<10) then
           do iapw=1,napwk(ikt)
              igv2(:,iapw,ikt)=igv2(:,iapw,ikt)+nlatout(:,1)
           enddo
           
        endif
      enddo
      deallocate(kv)
c      print *,'RRRRRRRRRRRRRRRRRRRRRRRRRRRRr suham.F'
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
c      if(iprx) 
      print *,' end of gen_hamindex: nkt,napwmx norbmto= ',nkt,napwmx,norbmto
      end subroutine gen_hamindex


      subroutine suham(sctrl,slat,sspec,ssite,sarray,sham,spot,sstrn,sbz)
      use m_struc_def  !Cgetarg
      use m_globalvariables
C- Hamiltonian setup
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Cio  sctrl :struct for program flow parameters; see routine uctrl
Cio    Elts read: nbasp nspec nl lham lcgf lgen3 lfp lncol nspin lpgf
Cio    Stored:    *
Cio    Passed to: *
Cio  slat  :struct for lattice information; see routine ulat
Cio    Elts read: ng ogv alat tolft plat qlat
Cio    Stored:    *
Cio    Passed to: *
Cio  sspec :struct for species-specific information; see routine uspec
Cio    Elts read: lmxa lmxl p pz
Cio    Stored:    *
Cio    Passed to: atfold makidx showbs sugcut
Cio  ssite :struct for site-specific information; see routine usite
Cio    Elts read: spec
Cio    Stored:    pnu pz
Cio    Passed to: showbs pvioeu
Cio  sarray:structure containing offsets to various arrays
Cio    Elts read: oips oipc nclasp opos
Cio    Stored:    *
Cio    Passed to: *
Cio  sham  :struct for parameters defining hamiltonian; see routine uham
Cio    Elts read: pwmode pwemax lncol oeula neula qss omagf nbf
Cio    Stored:    ndham ndofH ooffH oindxo ldham lmxax hord obdots
Cio    Passed to: *
Cio  spot  :struct for information about the potential; see routine upot
Cio    Elts read: opp vmtz
Cio    Stored:    nlma nlml opti
Cio    Passed to: *
Cio  sstrn :struct for global strings
Cio    Elts read: gfopt
Cio    Stored:    *
Cio    Passed to: *
Ci   gfopts:string containing switches for GF program
Cl Local variables
Cl   ndim  :total number of lmto orbitals = nl**2 * nbas
Cl   npwmin:lower limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   npwmax:upper limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   nqdiv: loop over mesh of q-points to estimate npwmax
Cl        : nqdiv is fineness of q-mesh.
Cl  npwpad: a 'safety' padding to npwmax in case npwmax
Cl        : underestimates actual upper limit
Cr Remarks
Cr   This routine generates energy-independent hamiltonian setup.
Cr  *It generates and packs a table of hamiltonian offsets offH,
Cr   orbital permutation indices oindxo.
Cr
Cr  *For the ASA 2nd generation LMTO:
Cr   Extract order of potential function from gfopts
Cr   Transform pp's to alpha representation
Cu Updates
Cu   07 Jul 08 Make sham->ndham = estimate for upper dimension of
Cu             hamiltonian, including possible PW part
Cu             Make sham->lmxax = largest lmax in basis
Cu   18 Apr 05 Force small parameter p -> 0 in 2-center turned on (ASA)
Cu   14 Feb 03 Makes and packs magnetic field B.sigma
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)
      type(s_array)::sarray
      type(s_ham)::sham
      type(s_pot)::spot
      type(s_bz)::sbz

      character*(*) sstrn
C ... Local parameters
      integer nkap0,n0H
      parameter (nkap0=3,n0H=5)
      integer hord,i,i1,i2,iprint,ldham(16),ldim, !lgen3,lham,
     .lgunit,lidim,lihdim,lncol,nbasp,nclasp,ndim,neul,nl,
     .nlspcp,nsp,nspc,nspec,nspx,nttab,partok,igets,nvi,nvl, !lfp,
     .ib,is,lmxa,lmxl,stdo,nglob,nkaph,isw,nbf,lmxax,pwmode,
     .j1,j2,j3,m,npw,npwmin,npwmax,ndham
      double precision pwemin,pwemax,plat(3,3),qlat(3,3),q(3),Gmin,Gmax,
     .xx
      integer n0,nqdiv,npwpad
      parameter (n0=10,nqdiv=12)
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))

C ... Needed for iostr
      double precision ckbas,cksumf,kap2(10)
      integer:: oalph , oiax , ontab , os , ng
ckino Dec.14.2011 manual correction       real(8),pointer :: rv_p_ogv(:) =>NULL()
      real(8),allocatable :: rv_a_ogv(:)

      integer ,allocatable :: kv_iv(:)
      integer ,allocatable :: igv2_iv(:)

      integer:: obas , oo , opp , oeula , opti , obs , omagf
ckino Dec.8.2011:            integer,pointer :: iv_p_oipc(:) =>NULL()

ckino Dec.8.2011:           integer,pointer :: iv_p_oips(:) =>NULL()

ckino Dec.14.2011:          integer,pointer :: iv_p_oidxsh(:) =>NULL()

ckino Dec.14.2011:         integer,pointer :: iv_p_offh(:) =>NULL()

      logical ltmp,iostr,bittst,adnf
      character*80 outs,gfopts
      double precision qss(4),vmtz,pnu(n0,2),pz(n0,2),alat,tolgv
C ... For NMTO
C     integer olmx,ohcr,otral,otrad,oalpha,oadot,nmto
C     double precision kmto(10),avw,dglob
C ... Heap

c      integer :: sigswitch

C --- Setup ---

ckino Dec.8.2011:         iv_p_oips => sarray%iv_p_oips

ckino Dec.8.2011:         iv_p_oipc => sarray%iv_p_oipc



      nbasp=sctrl%nbasp
      nspec=sctrl%nspec
      nl=sctrl%nl


c       lham=sctrl%lham
c       lcgf=sctrl%lcgf
c       lgen3=sctrl%lgen3
c       lfp=sctrl%lfp


      lncol=sctrl%lncol
      nsp=sctrl%nspin

c       lpgf = int(sctrl%lpgf(1))

Changenglob      nkaph = nglob('nkaph')
      nkaph = globalvariables%nkaph

      pwmode=sham%pwmode
      pwemin=sham%pwemin
      pwemax=sham%pwemax
      npwpad=sham%npwpad


      nspc  = 1
      nspx  = nsp
      if (lncol .ne. 0) then
        nspc = 2
        nspx = 1
      endif
      ndim = nbasp * nl**2 * nkaph
ckino Dec.14.2011:         allocate(iv_p_offh(abs(-n0h*nkap0*(nbasp+1))))
ckino Dec.28.2011:         allocate(sham%iv_p_ooffh(abs(-n0h*nkap0*(nbasp+1))))
      if (allocated(sham%iv_a_ooffH)) deallocate(sham%iv_a_ooffH)
      allocate(sham%iv_a_ooffH(abs(-n0h*nkap0*(nbasp+1))))
ckino Dec.14.2011:         if (-n0h*nkap0*(nbasp+1)<0) iv_p_offh(:)=0
ckino Dec.28.2011:         if (-n0h*nkap0*(nbasp+1)<0) sham%iv_p_ooffh(:)=0
      if (-n0h*nkap0*(nbasp+1)<0) sham%iv_a_ooffH(:)=0

ckino Dec.14.2011:         allocate(iv_p_oidxsh(abs(ndim+3)))
ckino Dec.27.2011:         allocate(sham%iv_p_oindxo(abs(ndim+3)))
      allocate(sham%iv_a_oindxo(abs(ndim+3)))
ckino Dec.14.2011:         if (ndim+3<0) iv_p_oidxsh(:)=0
ckino Dec.27.2011:         if (ndim+3<0) sham%iv_p_oindxo(:)=0
      if (ndim+3<0) sham%iv_a_oindxo(:)=0

      stdo = lgunit(1)

C --- Automatic downfolding (2nd generation only?) ---
c      if (lgen3 .eq. 0 .and. lfp .eq. 0) then
C$$$      if (lfp .eq. 0) then
C$$$        adnf = bittst(lham,4)
C$$$C       Suppress if in gamma repsn (not needed for GF mode)
C$$$C       if (lgors('ctrl lham,128',sctrl)) adnf = 0
C$$$Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
C$$$Cgetarg...info...           real(8) :: vmtz , real(8) :: spot%vmtz 1
C$$$Cgetarg         call upack2('pot opp vmtz',spot,opp,vmtz)

C$$$         opp=spot%opp
C$$$         vmtz=spot%vmtz

C$$$        call atfold(0,adnf,nl,nsp,vmtz,w(opp),sspec)
C$$$      endif

C --- Hamiltonian offsets, orbital permutation table ---
      call iinit(ldham,16)
      if (mod(pwmode,10) .eq. 2) then
        call info0(20,0,0,' suham: LMTO basis will be excluded')
      else
ckino Dec.8.2011:             call makidx ( nl , nkaph , 1 , nbasp , 0 , sspec , iv_p_oips
ckino Dec.28.2011:             call makidx ( nl , nkaph , 1 , nbasp , 0 , sspec , sarray%iv_p_oips
        call makidx ( nl , nkaph , 1 , nbasp , 0 , sspec , sarray%iv_a_oips
ckino Dec.14.2011:       .     , iv_p_offh , iv_p_oidxsh , ldham )
ckino Dec.14.2011:       .     , iv_p_offh , sham%iv_p_oindxo , ldham )
ckino Dec.27.2011:       .     , sham%iv_p_ooffh , sham%iv_p_oindxo , ldham )
ckino Dec.28.2011:       .     , sham%iv_p_ooffh , sham%iv_a_oindxo , ldham ) 
     .     , sham%iv_a_ooffH , sham%iv_a_oindxo , ldham )



      endif
C     Default value for hamiltonian dimension

      sham%ndham=ldham(1)

C     print *, '!!'; call pack1('ham ndham',sham,ldham+3)

C ... Make rest of ldham
      ldham(4) = nspc
      ldham(5) = ldham(1) * nspc
      ldham(6) = ldham(2) * nspc
      ldham(7) = ldham(3) * nspc
      ldham(8) = nspx

      sham%ndofH=nbasp+1
ckino Dec.14.2011:         sham%iv_p_ooffh => iv_p_offh

ckino Dec.14.2011:         sham%iv_p_oindxo => iv_p_oidxsh

      i_copy_size=size(sham%ldham)
      call ii8copy(i_copy_size,ldham,1,sham%ldham,1)


C ------- Potential- and implementation- specific setup -------
C --- FP setup ---

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)

c      if (lfp .ne. 0) then
      nvi = 0
      nvl = 0
      lmxax = -1
      do  20  ib = 1, nbasp

        is=ssite(ib)%spec


        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl
        i_copy_size=size(sspec(is)%p)
        call dcopy(i_copy_size,sspec(is)%p,1,pnu,1)
        i_copy_size=size(sspec(is)%pz)
        call dcopy(i_copy_size,sspec(is)%pz,1,pz,1)

cccccccccccccc
c        print *,'suham 111 pnu=',i,pnu(1:10,1:2)
cccccccccccccc
        call dcopy(n0,pz,1,pz(1,2),1)
C         Augmentation dimensioning parameters
        nvi = nvi + (lmxa+1)**2
        nvl = nvl + (lmxl+1)**2
C         Poke spec starting pnu to site

        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,pz,1,ssite(ib)%pz,1)

C         Find largest lmxa
        lmxax = max(lmxax,lmxa)
   20 continue

      sham%lmxax=lmxax


      spot%nlma=nvi
      spot%nlml=nvl

      call info5(30,0,0,' suham :  %i augmentation'//
     .' channels, %i local potential channels  Maximum lmxa=%i',
     .nvi,nvl,lmxax,0,0)

      ng=slat%ng
ckino Dec.14.2011       rv_p_ogv => slat%rv_p_ogv

      alat=slat%alat
      tolgv=slat%tolft

ckino Dec.14.2011       call sugcut ( 1 , nspec , sspec , alat , ng , rv_p_ogv , tolgv
ckino Dec.27.2011:         call sugcut ( 1 , nspec , sspec , alat , ng , slat%rv_p_ogv , tolgv
      call sugcut ( 1 , nspec , sspec , alat , ng , slat%rv_a_ogv , tolgv 
     . )


C   ... PW setup : estimate upper bound to number of G vectors
C       to set up upper bound to hamiltonian dimension

      pwmode=sham%pwmode
      pwemax=sham%pwemax

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
          if (npwpad .lt. 0) then
            npwpad = max(nint((npwmax-npwmin)*0.2d0),3)
          endif
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

        sham%npwmin=npwmin
        sham%npwpad=npwpad
        sham%ndham=ndham

        if (mod(pwmode/10,10) .eq. 1) then
          call info2(20,1,0,' suham:  q-dependent PW basis with'//
     .    '  Emin = %d < E < %d.',pwemin,pwemax)
          call info5(30,0,0,'%9fEst. min,max PW dimension = %i,%i.'//
     .    '  Use npwpad = %i => ndham = %i',
     .    npwmin,npwmax,npwpad,ndham,0)
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

ckino Dec.14.2011 manual correction           allocate(rv_p_ogv(abs(3*npw)))
ckino Dec.14.2011 manual correction           if (3*npw<0) rv_p_ogv(:)=0.0d0
          allocate(rv_a_ogv(abs(3*npw))); rv_a_ogv(:)=0.0d0

          call pshpr(iprint())
          if (iprint() .ge. 50) call setpr(100)
          call gvlst2 ( alat , plat , q , 0 , 0 , 0 , gmin , gmax , 0 ,
ckino Dec.14.2011 manual correction     .     8 + 2 , npw , npw , kv_iv , rv_p_ogv , xx , igv2_iv ) 
     .     8 + 2 , npw , npw , kv_iv , rv_a_ogv , xx , igv2_iv )


          call poppr
          if (allocated(kv_iv)) deallocate(kv_iv)
          if (allocated(igv2_iv)) deallocate(igv2_iv)
ckino Dec.14.2011 manual correction           if (associated(rv_p_ogv)) deallocate(rv_p_ogv)
          deallocate(rv_a_ogv)


        endif
      endif
c      endif

C ... Printout orbital positions in hamiltonian, resolved by l
c      if (iprint() .gt. 50) then
c        call showbs(ssite,sspec,slat,nkaph,nl,0,w(oidxsh),ldham)
c      endif




c this must be used for the case of QSGW.
c      if(sigswitch()>0) then
c        call gen_hamindex(sham,sbz,ssite,sspec,slat,nl,nbasp,w(oidxsh),ldim)
c      endif


C --- Third-generation NMTO ---
c      if (lgen3 .ne. 0) then
c         stop 'xxxxxxxxxxx third generation nmto path. not meaningful'
C        call upack('ham nmto kmto',sham,nmto,kmto,0,0,0)
Ckino isanrg is logical function, C        call isanrg(nmto,2,5,'suham','nmto',.true.)
C        l_dummy_isanrg=isanrg(nmto,2,5,'suham','nmto',.true.)
C        call defdr(otral,-4*nl**2*nbasp*nmto)
C        call defdr(otrad,-4*nl**2*nbasp*nmto)
C
C        call sp2cls('spec lmxb',sspec,0,1,1,nspec,olmx)
C        call sp2cls('spec hcr',sspec,0,nl,1,nspec,ohcr)
C        call defdr(oalpha,-nl**2*nbasp*nmto)
C        call defdr(oadot,-nl**2*nbasp*nmto)
C        avw = dglob('avw',0d0,0)
C        call dscal(nmto,avw**2,kmto,1)
C        call dscal(nl*nspec,1/avw,w(ohcr),1)
C        call pshpr(50)
C        call mktra2(1,loka,nbasp,w(oips),nl,w(olmx),avw,4,kmto,nmto,
C     .    w(ohcr),1,w(otral),w(otrad),w(oalpha),w(oadot),w,w)
C        call poppr
C
C        call dscal(nmto,1/avw**2,kmto,1)
C        call dscal(nl*nspec,avw,w(ohcr),1)
C        call rlse(olmx)
c      endif

C --- Second-generation LMTO, ASA ---
C ... Green's function-specific initialization
C$$$      if (lcgf .ne. 0 .or. lpgf .ne. 0) then
C$$$        call upacks('strn gfopt',i1,i2)
C$$$        gfopts = ' '
C$$$        if (i2 .ge. i1) then
C$$$          gfopts = sstrn(i1:i2)
C$$$          call partk0(0,len(gfopts),1,-1,0,len(gfopts),-1,31,.false.)
C$$$          hord = 2
C$$$          ltmp = .false.
C$$$          i = partok(gfopts,'pz',  ' ;',ltmp,' ',0,0,0,0)
C$$$          if (ltmp) hord = 4
C$$$          i = partok(gfopts,'p3',  ' ;',ltmp,' ',0,0,0,0)
C$$$          if (ltmp) hord = 3
C$$$          i = partok(gfopts,'p1',  ' ;',ltmp,' ',0,0,0,0)
C$$$          if (ltmp) hord = 1
C$$$Cgetarg...info...           integer :: hord , integer(8) :: sham%hord 1
C$$$Cgetarg           call pack1('ham hord',sham,hord)
C$$$           sham%hord=hord
C$$$          if (iprint() .ge. 30) call awrit1(' SUHAM : 2nd generation'//
C$$$     .      ' ASA; potential functions P%?#n<4#%-1j%i#(z)#',' ',
C$$$     .      80,stdo,hord)
C$$$        endif
C$$$      endif

C ... 2nd gen LMTO-specific initialization
C     Make pti here since pp in gamma rep
c      if (lcgf .eq. 0 .and. lfp .eq. 0 .and. lgen3 .eq. 0) then
c      if (lfp .eq. 0 .and. lgen3 .eq. 0) then
C$$$      if (lfp .eq. 0) then

C$$$        call defdr(opti,-ndim*nsp)
C$$$        call makipt(nl,nbasp,nsp,w(oipc),lihdim,w(oidxsh),w(opp),ndim,
C$$$     .    w(opti))
C$$$Cgetarg...info...           integer :: opti , integer(8) :: spot%opti 1
C$$$Cgetarg         call pack1('pot opti',spot,opti)

C$$$         spot%opti=opti


C$$$      endif

C ... Transform pp's to alpha representation
c      if (lfp .eq. 0 .and. lgen3 .eq. 0) then
C$$$      if (lfp .eq. 0) then
C$$$C       Get screening alpha's
C$$$Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
C$$$Cgetarg...info...           integer :: obas , integer(8) :: sarray%opos 1
C$$$Cgetarg         call upack('array nclasp opos',sarray,nclasp,obas,0,0,0)

C$$$         nclasp=sarray%nclasp
C$$$         obas=sarray%opos

C$$$        nlspcp = nl * nsp * nclasp
C$$$        ckbas = cksumf(w(obas),3*nbasp)
C$$$        ltmp = iostr(8,'STR',nl,nbasp,1,kap2,0,ckbas,-1,nttab,oalph,
C$$$     .    oiax,ontab,os)
C$$$        i = 0
C$$$        if (bittst(lham,1) .or. bittst(lham,2)) i = 10
C$$$C       Do the transformation
C$$$        call defdr(oo,nlspcp)
C$$$        call pptrns(i,nl,w(oipc),nclasp,nsp,w(oalph),1,w(opp),w(oo))
C$$$        call rlse(oalph)
C$$$      endif

C ... Printout Euler angles

c       lncol=sham%lncol
c       oeula=sham%oeula
c       neul=sham%neula
c       i_copy_size=size(sham%qss)
c       call dcopy(i_copy_size,sham%qss,1,qss,1)

C$$$      if (iprint() .ge. 30 .and. lncol .ne. 0) then
C$$$        call awrit4('%x          %?#n#Qss%3:1;6d angle %;6d.  #%2j#'//
C$$$     .    '%?#n#Euler angles:#noncollinear hamiltonian',outs,80,
C$$$     .    -stdo,bittst(lncol,2),qss,qss(4),isw(bittst(lncol,1).and.
C$$$     .    iprint().gt.30))
C$$$        if (iprint().gt.30 .and. bittst(lncol,1))
C$$$     .    call pvioeu(1,ssite,w(oeula),nbasp,neul)
C$$$      endif

C ... Make B.sigma for magnetic field
c      obs = 1
c      if (lfp.eq.0 .and. lgen3.eq.0 .and. bittst(lncol,8)) then
C$$$      if (lfp.eq.0 .and. bittst(lncol,8)) then
C$$$Cgetarg...info...           integer :: oeula , integer(8) :: sham%oeula 1
C$$$Cgetarg...info...           integer :: omagf , integer(8) :: sham%omagf 1
C$$$Cgetarg...info...           integer :: nbf , integer(8) :: sham%nbf 1
C$$$Cgetarg         call upack('ham oeula omagf nbf',sham,oeula,omagf,nbf,0,0)

C$$$         oeula=sham%oeula
C$$$         omagf=sham%omagf
C$$$         nbf=sham%nbf

C$$$        call defcc(obs,-4*lihdim)
C$$$        call mkbfld(nl,nbasp,lihdim,w(oidxsh),w(omagf),nbf,w(oeula),
C$$$     .    neul,w(obs))
C$$$      endif

c       sham%obdots=obs


      end subroutine suham






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
      integer  iprint
      external iprint
      integer ires(3, nbas, ng)
ckino delete integer(4) def.      integer(4):: ib1,ib2
      integer:: ib1,ib2
c
      real(8) ::tran(3),delta(3,ng)
      data ep/1.0d-3/
c      data ep/1.0d-7/
c
      if(iprint()>=46) write(6,*)'MPTAUOf: search miat tiat for wave function rotation'

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
        if(iprint() .ge.46) then
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
        stop 'mptauof: Can not find miat and tiat'
c
   21   continue
        delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000

c- have gotten the translation-> check write --------------------
        if(iprint().ge.46) then
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

      igann  = 1d0
      do ix =1,i
        igann=igann*dble(ix)
      enddo
      end function igann

