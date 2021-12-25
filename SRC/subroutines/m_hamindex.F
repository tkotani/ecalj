      module m_hamindex
      use NaNum,only: NaN       !not working well
      integer,protected:: ngrpaf,ngrp_original
      integer,allocatable,protected:: iclasstaf(:)
      integer, allocatable,protected ::  offH (:,:,:)
      real(8),allocatable,protected:: symops_af(:,:,:), ag_af(:,:)
! we have to protect them in future      
      integer,allocatable:: !,protected::
     &     ltab(:),ktab(:),offl(:),ispec(:), iclasst(:),offlrev(:,:,:),ibastab(:),
     &     iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:), !,ngvecp(:,:,:),ngvecprev(:,:,:,:)
     &     igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
      real(8):: plat(3,3)=NaN,qlat(3,3)=NaN !,protected
      real(8),allocatable:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:),
     &     qtt(:,:),qtti(:,:)   !,protected
      integer::  nqi=NaN, nqnum=NaN, ngrp=NaN, lxx=NaN, kxx=NaN,norbmto=NaN,
     &      nqtt=NaN, ndimham=NaN, napwmx=NaN, lxxa=NaN, ngpmx=NaN !, norbtx=NaN
      integer:: imx=NaN,nbas=NaN !,protected
      logical,private:: debug=.false.
      contains
cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine m_hamindex_init(ljobgw)
      use m_mksym,only: rv_a_osymgr,rv_a_oag,lat_nsgrp, iclasstaf_,symops_af_,ag_af_,ngrpaf_
      use m_struc_def           !Cgetarg
      use m_lmfinit,only:
     &     ham_pwmode,pwemax,ham_ldham,
     &     alat=>lat_alat,nl,ctrl_nbas,ssite=>v_ssite,sspec=>v_sspec,n0,nkap0,iprmb
      use m_lattic,only: lat_qlat,lat_plat,rv_a_opos
      use m_shortn3,only: shortn3_initialize,shortn3
      use m_MPItk,only: master_mpi
!!-- Set up m_hamiltonian. Index for Hamiltonian. --
!!  Generated index are stored into m_hamindex
!!  Only include q-point information for GW (QGpsi).
!!#Inputs
Cr As you see in subroutine rotwv, the index for Hamiltonian reads as;
c      do iorb=1,norbmto             !orbital-blocks are specified by ibas, l, and k.
c        ibas  = ibastab(iorb)
c         l    = ltab(iorb)
c         k    = ktab(iorb)        !kappa (diffent MTO index for each l)
c        init1 = offl(iorb)+1      !starting index for the block iorb
c        iend1 = offl(iorb)+2*l+1  !end of the block for the iorb
c      enddo
      implicit none
      integer,parameter:: noutmx=48, n0nkap0=30, n0H=5
      integer:: ibas,k,l,ndim,ipr,nglob,ldim,off,offs,specw,fieldw,iorb,offsi,ib,is, norb
      integer:: nkabc(3),nkp,lshft(3),pwmode,napwx,ig,nini,nk1,nk2,nk3,ik1,ik2,ik3,ikt,nkt
      integer:: i_copy_size,i_spacks,i_spackv,ifi,nbas_in,ifisym,i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw
      integer:: iout,nout,nlatout(3,noutmx),iapw ,iprint,ngadd,igadd,igaf
      integer:: ngp, ifiqg,iq,fopnx,nnn(3),ixx,ndummy,nqbz___ ,ifatomlist
      integer,allocatable:: ltabx(:,:),ktabx(:,:),offlx(:,:),iqtt(:), kv(:)
      real(8):: pwgmax, pwgmin, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)
      real(8):: dum,qb(3,3),ddd(3),ppin(3), tolq, rlatp(3,3),xmx2(3),qqq(3),diffs,ddf
      real(8),allocatable:: symtmp(:,:,:)
      logical:: siginit, qpgexist,debug=.false., llmfgw,ljobgw, prpushed
      character(8),allocatable::  spid(:)
!!      
      ldim  = ham_ldham(1)
      nbas=ctrl_nbas !nbas_in
      ngrp=lat_nsgrp !note nsgrp given in mksym.F is without inversion.
      plat=lat_plat
      qlat=lat_qlat
      allocate( ltabx(n0nkap0,nbas),ktabx(n0nkap0,nbas),offlx(n0nkap0,nbas)) !,ispec(nbas))
c      if(master_mpi) write(6,*) ' gen_hamindex: --- Hamiltonian index ---'
c      if(master_mpi) write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
!!--- MTO part ---
!!---- obtain norbmto, lxx,kxx ----
      norbmto=0
      kxx=-1
      lxx=-1
      ndimham = 0 !dimension of mto part of hamiltonian
      do  ib = 1, nbas
        is=ssite(ib)%spec
        call orbl(ib,0,ldim,iprmb,norb,ltabx(:,ib),ktabx(:,ib),off,offlx(:,ib),ndim)!iprmb
        do iorb = 1, norb
          norbmto = norbmto+1
          if(ltabx(iorb,ib)>lxx)  lxx = ltabx(iorb,ib)
          if(ktabx(iorb,ib)>kxx)  kxx = ktabx(iorb,ib)
          ndimham = ndimham+ 2*ltabx(iorb,ib)+1
        enddo
      enddo
!!--- make index table :norbmto is the total number of different type of MTOs 
      allocate( ibasindex(ndimham))
      allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
      allocate( spid(nbas))
      norbmto=0
      ndimham = 0 !dimension of mto part of hamiltonian
      allocate(offH(n0h,nkap0,nbas+1))
      offH=0
      do  ib = 1, nbas
        is=ssite(ib)%spec
        spid(ib)=sspec(is)%name
        call orbl(ib,0,ldim,iprmb,norb,ltabx(:,ib),ktabx(:,ib),off,offlx(:,ib),ndim) !iprmb
        do  iorb = 1, norb
          norbmto=norbmto+1
          ibastab(norbmto)= ib
          ltab(norbmto)   = ltabx(iorb,ib) !angular momentum l index
          ktab(norbmto)   = ktabx(iorb,ib) !
          offl(norbmto)   = offlx(iorb,ib)
          nini = ndimham+ 1
          ndimham = ndimham+ 2*ltabx(iorb,ib)+1
          ibasindex(nini:ndimham) = ib
c          write(6,"('offx ',3i3,2x,2i5,3x,a)") 
c     &    ib,ltab(norbmto),ktab(norbmto), offl(norbmto)+1,ndimham,trim(spid(ib))
        enddo
        offH(1,1,ib+1) = ndimham
      enddo
      offH(1,1,nbas+1) = ndimham
!!      
      ifatomlist=999
      open(ifatomlist,file="ATOM_LIST")
      write(ifatomlist,*)nbas
      write(ifatomlist,*)spid
      close(ifatomlist)
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
      if (master_mpi) then
      open(newunit=ifisym,file='SYMOPS')
      write(ifisym,*) ngrp
      do ig = 1,ngrp
        write(ifisym,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
        do i=1,3
          write(ifisym,"(3e24.16)") symops(i,1:3,ig)
        enddo
      enddo
      close(ifisym)
      endif
!! Get rotation matrix dlmm.  We assume nl=lmxa+1.
      lxxa=nl-1
      allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      call rotdlmm(symops,ngrp, nl, dlmm)

!! No GW-related part
      if(ljobgw) then
c      if((.not.llmfgw).or.(llmfgw.and.jobgw==0)) then
        if(master_mpi) print *,'gen_hamindex: not readin QGpsi.'
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

!! we have two set of data for original qxx in QGpsi and their shortened.
c      nqnum2 = 2*nqnum !+ 2*nkp !takao feb2012 test xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(allocated(qq)) deallocate(qq)
      nqtt=nqnum !nqnum*2 !doubled. second series nqnum+1:2*nqnum are for shortened q.
      nkt=nqtt
      allocate( qtti(3,nqi), qq(3,nqtt),iqtt(nqtt) )
      if(master_mpi) print *
      if(master_mpi) print *,'gen_hamindex: Readin QGpsi. nqnum=',nqnum
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
        if(master_mpi)write(6,"(' qq=',i5,3f10.5)") iq,qq(:,iq)
      enddo
      close(ifiqg)

!! ==== Generate info for rotwv and write ====
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

        if(master_mpi) then
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
c      pwemax=ham_pwemax
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
      if(master_mpi) print *,' gen_hamindex goto APW part: pwmode pwemax=',pwmode,pwemax !pwemin
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
      if(master_mpi)print*,' --- gvlst2 generates G for APW part (we show cases for limited q) ---'
      if(pwmode<5) call shortn3_initialize(qlat)
      allocate( igv2(3,napwmx,nkt), kv(3*napwmx) )
      prpushed=.false.
      do ikt = 1,nkt
        qqq=0d0
        if (mod(pwmode/10,10) == 1) qqq=qq(:,ikt) 
        call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napwmx,napwk(ikt),kv,dum,dum,igv2(1,1,ikt))
        if(master_mpi.and.(ikt>5.or.ikt==nkt).and.(.not.prpushed)) then
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
      if(master_mpi) call poppr !print index is poped.
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
      if(master_mpi) print *,' end of gen_hamindex: nkt,napwmx norbmto= ',nkt,napwmx,norbmto
!! Symmetry for AF. Order AF symmetry operation after normal one. jun2015takao
!! Caution: symops,and so on are overwritten.
!!          writehamindex() already wrote HAMindex which is just for SYMGRP.
 2001 continue
      if(allocated(symops_af)) then
        ngrpaf = ngrpaf_
        allocate(iclasstaf(nbas),symops_af(3,3,ngrpaf_),ag_af(3,ngrpaf_))
        iclasstaf = iclasstaf_
        symops_af = symops_af_
        ag_af = ag_af_
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
        symtmp(:,:,igadd)=symops_af(:,:,igaf)
 1013   continue
        enddo
        if(igadd/=ngrpaf) call rx('suham: strange. bug igadd/=ngrpaf')
        if(master_mpi) write(6,*) '-----SYMGRPAF mode ---- # of additional symmetry=',igadd
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
        if(master_mpi) then
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
      endif
      end subroutine m_hamindex_init
      
      
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
!> get index ikt such that for qin(:)=qq(:,ikt)
      integer function getikt(qin) !return
      intent(in)::            qin
      integer::i
      real(8):: qin(3)
c      if(debug) print *,'nkt=',nkt
      do i=1, nqnum !*2 !nkt
        if(debug) print *,i,qin, qq(:,i)
        if(sum (abs(qin-qq(:,i)))<1d-8) then
          getikt=i
          return
        endif
      enddo
      print *,' getikt: xxx error nqnum qin=',nqnum,qin
      do i=1, nqnum !*2 !nkt
        write(*,"('i qq=',i3,3f11.5)")i, qq(:,i)
      enddo  
      call rx( ' getikt can not find ikt for given q')
      end function

!> write info for wave rotation. 
      subroutine writehamindex()
      integer(4):: ifi
      logical::pmton
      logical,save:: done=.false.
      if(done) call rx('writehamindex is already done')
      done=.true.
      open(newunit=ifi,file='HAMindex',form='unformatted')
      write(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto
      write(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
      write(ifi)lxxa
      write(ifi)dlmm
      write(ifi)ibastab,ltab,ktab,offl,offlrev !for rotation of MTO. recovered sep2012 for EIBZ for hsfp0
      write(ifi)qq !,ngvecp,ngvecprev
      write(ifi)plat,qlat,napwmx
      if(napwmx/=0) then !for APW rotation used in rotwvigg
         write(ifi) igv2,napwk,igv2rev
      endif   
      write(ifi)plat,qlat,napwmx
      close(ifi)
      end subroutine writehamindex

!> read info for wave rotation. 
      subroutine readhamindex()
!! == read info for wave rotation. feb2012takao ==
      integer(4):: ifi,nkt
      logical::pmton
      logical,save:: done=.false.
      if(done) call rx('readhamindex is already done')
      done=.true.
      open(newunit=ifi,file='HAMindex',form='unformatted')
      read(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto
      allocate(symops(3,3,ngrp),ag(3,ngrp),qtt(3,nqtt),qtti(3,nqi))
      allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
      write(6,*) 'ngrp=',ngrp
      read(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
      allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
      allocate( offlrev(nbas,0:lxx,kxx))
      read(ifi) lxxa
      allocate( dlmm(-lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      read(ifi) dlmm
      read(ifi)ibastab,ltab,ktab,offl,offlrev
c      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
c      allocate( ngvecp(3,ngpmx,nqnum) )
      allocate( qq(3,nqnum)) !this was qq(3,nqnum*2) until Aug2012 when shorbz had been used.
      read(ifi)qq !,ngvecp,ngvecprev
      read(ifi)plat,qlat,napwmx
      if(napwmx/=0)then !for APW rotation used in rotwvigg
        nkt=nqtt
        allocate( igv2(3,napwmx,nkt) )
        allocate( napwk(nkt))
        allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt) )
        read(ifi) igv2,napwk,igv2rev
      endif  
      close(ifi)
      end subroutine readhamindex
      
      end module


