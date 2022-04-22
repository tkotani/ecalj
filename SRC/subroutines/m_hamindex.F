      module m_hamindex         !all protected now 
      use m_lmfinit,only: ham_pwmode,pwemax,ham_ldham,noutmx,nsp_in=>nsp,stdo,
     &     alat=>lat_alat,nl,ctrl_nbas,ssite=>v_ssite,sspec=>v_sspec,n0,nkap0,iprmb,zbak_read=>zbak
      use m_lattic,only: lat_qlat,lat_plat,rv_a_opos
      use NaNum,only: NaN       !for initialization, but not working well
      use m_suham,only: ndham_read=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
      integer,protected,public:: ngrpaf,ngrp_original,pwmode,ndham
      integer,protected,public:: nqi=NaN, nqnum=NaN, ngrp=NaN, lxx=NaN, kxx=NaN,norbmto=NaN,
     &     nqtt=NaN, ndimham=NaN, napwmx=NaN, lxxa=NaN, ngpmx=NaN, imx=NaN,nbas=NaN
      integer,allocatable,protected,public:: iclasstaf(:), offH (:),
     &     ltab(:),ktab(:),offl(:),ispec(:), iclasst(:),offlrev(:,:,:),ibastab(:),
     &     iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:),
     &     igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
      real(8),allocatable,protected,public:: symops_af(:,:,:), ag_af(:,:),
     &     symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:),
     &     qtt(:,:),qtti(:,:)  
      real(8),protected,public:: plat(3,3)=NaN,qlat(3,3)=NaN,zbak 
      public:: m_hamindex_init, Readhamindex, getikt

      private
      logical,private:: debug=.false.
      contains
cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine m_hamindex_init(jobgw)
      use m_mksym,only: rv_a_osymgr,rv_a_oag,lat_nsgrp, iclasstaf_,symops_af_,ag_af_,ngrpaf_
      use m_struc_def           !Cgetarg
      use m_shortn3,only: shortn3_initialize,shortn3
      use m_MPItk,only: master_mpi
!!-- Set up m_hamiltonian. Index for Hamiltonian. --
!!  Generated index are stored into m_hamindex
!!  Only include q-point information for GW (QGpsi).
!!#Inputs
Cr As you see in subroutine rotwvigg, the index for Hamiltonian reads as;
c      do iorb=1,norbmto             !orbital-blocks are specified by ibas, l, and k.
c        ibas  = ibastab(iorb)
c         l    = ltab(iorb)
c         k    = ktab(iorb)        !kappa (diffent MTO index for each l)
c        init1 = offl(iorb)+1      !starting index for the block iorb
c        iend1 = offl(iorb)+2*l+1  !end of the block for the iorb
c      enddo
      implicit none
      integer:: ibas,k,l,ndim,ipr,nglob,ldim,off,offs,specw,fieldw,iorb,offsi,ib,is, norb,nsp
      integer:: nkabc(3),nkp,lshft(3),napwx,ig,nini,nk1,nk2,nk3,ik1,ik2,ik3,ikt,nkt
      integer:: i_copy_size,i_spacks,i_spackv,ifi,nbas_in,ifisym,i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw
      integer:: iout,nout,nlatout(3,noutmx),iapw ,iprint,ngadd,igadd,igaf
      integer:: ngp, ifiqg,iq,nnn(3),ixx,ndummy,nqbz___ ,ifatomlist
      integer,allocatable:: ltabx(:,:),ktabx(:,:),offlx(:,:),iqtt(:), kv(:)
      real(8):: pwgmax, pwgmin, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)
      real(8):: dum,qb(3,3),ddd(3),ppin(3), tolq, rlatp(3,3),xmx2(3),qqq(3),diffs,ddf
      real(8),allocatable:: symtmp(:,:,:)
      logical:: siginit, qpgexist,debug=.false., llmfgw,prpushed
      character(8)::  spid(ctrl_nbas)
      integer::ndima,lmxax,npqn,ificlass,nat,lmaxa,ipqn,ifinlaindx,isp,konf
      integer,allocatable:: konft(:,:,:),iqnum(:,:)
      real(8) ::xx(5),pnu(n0,2),pnz(n0,2)
      integer:: ipb(ctrl_nbas),ipc(ctrl_nbas),ipcx(ctrl_nbas),mxint
      character lsym(0:n0-1)*1, lorb(3)*1, dig(9)*1, strn4*4
      character(4),allocatable:: strn4c(:)
      data dig /'1','2','3','4','5','6','7','8','9'/
      data lsym /'s','p','d','f','g','5','6','7','8','9'/
      data lorb /'p','d','l'/
! forget floating orbital case
      ldim  = ham_ldham(1) 
      nbas=ctrl_nbas !nbas_in
      ngrp=lat_nsgrp !note nsgrp given in mksym.F is without inversion.
      plat=lat_plat
      qlat=lat_qlat
      nsp=nsp_in
!! symmetry operation ---
      allocate(symops(3,3,ngrp),ag(3,ngrp))
      call dcopy ( ngrp * 9 , rv_a_osymgr , 1 , symops , 1 )
      call dcopy ( ngrp * 3 , rv_a_oag , 1 , ag , 1 )
      allocate(iclasst(nbas),invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      do ib=1,nbas
        iclasst(ib)=ssite(ib)%class
      enddo
!! get space group information ---- translation informations also in miat tiat invgx, shtvg
      call mptauof( symops , ngrp , plat , nbas , rv_a_opos, iclasst 
     .   , miat , tiat , invgx , shtvg )
c$$$!! jobgw=0 mode ----------------------------------
c$$$      if(master_mpi.and.jobgw==0) then
c$$$         ndima = 0
c$$$         do  ib = 1, nbas
c$$$           is=ssite(ib)%spec
c$$$           lmxa(ib) = sspec(is)%lmxa  !we assume lmxa>-1 
c$$$           call dcopy(size(sspec(is)%pz),sspec(is)%pz,1,pnz,1)
c$$$c           if (lmxa(ib)> -1) then
c$$$              do  l = 0, lmxa(ib)
c$$$                npqn = 2
c$$$                if (pnz(l+1,1) .ne. 0) npqn = 3
c$$$                ndima = ndima + npqn*(2*l+1)
c$$$              enddo
c$$$c           endif
c$$$         enddo
c$$$         lmxax = mxint(nbas,lmxa)
c$$$         allocate(konft(0:lmxax,nbas,nsp))
c$$$         do ib = 1, nbas
c$$$           call dcopy(size(ssite(ib)%pnu), ssite(ib)%pnu,1,pnu,1)
c$$$           call dcopy(size(ssite(ib)%pz),  ssite(ib)%pz, 1,pnz,1)
c$$$           do  isp = 1, nsp
c$$$             do  l  = 0, lmxa(ib)
c$$$               konft(l,ib,isp) = pnu(l+1,isp)
c$$$               if( mod(pnz(l+1,isp),10d0)<pnu(l+1,isp) .and. pnz(l+1,isp)>0) then
c$$$                  konft(l,ib,isp) = mod(pnz(l+1,isp),10d0)
c$$$               endif    
c$$$             enddo
c$$$           enddo
c$$$         enddo
c$$$!     ! NLAindx
c$$$c         open(newunit=ifinlaindx,file='NLAindx')
c$$$c         write(ifinlaindx,'(''----NLAindx start---------------''/I6)') ndima
c$$$         ndima = 0
c$$$         allocate(iqnum(4,3*nbas*lmxax),strn4c(3*nbas*lmxax))
c$$$         iorb=0
c$$$         do  ipqn = 1, 3
c$$$         do  ib = 1, nbas
c$$$             is =  ssite(ib)%spec
c$$$             lmaxa=sspec(is)%lmxa
c$$$             i_copy_size=size(ssite(ib)%pnu)
c$$$             call dcopy(i_copy_size ,ssite(ib)%pnu,1,pnu,1)
c$$$             i_copy_size=size(ssite(ib)%pnu)
c$$$             call dcopy(i_copy_size,ssite(ib)%pz, 1,pnz,1)
c$$$c             if (lmaxa .gt. -1) then
c$$$                do  l = 0, lmaxa
c$$$                  npqn = 2
c$$$                  if (pnz(l+1,1) .ne. 0) npqn = 3
c$$$                  if (ipqn <= npqn) then
c$$$                     iorb=iorb+1
c$$$                     konf = pnu(l+1,1)
c$$$                     if (ipqn .eq. 3) konf = mod(pnz(l+1,1),10d0)
c$$$                     strn4 = dig(konf)//lsym(l)//'_'//lorb(ipqn)
c$$$                     iqnum(:,iorb)=[ipqn,l,ib,ndima]
c$$$                     strn4c(iorb)=strn4
c$$$c                     write(ifinlaindx,'(i6,i3,i4,i6,4x,a)')ipqn,l,ib,ndima,strn4
c$$$c                    write(ifinlaindx,'(i6,i3,i4,i6,4x,a)')ipqn,l,ipb(ib),ndima,strn4
c$$$                     ndima = ndima + (2*l+1)
c$$$                  endif
c$$$                enddo
c$$$c             endif
c$$$         enddo
c$$$         enddo
c$$$         norb=iorb
c$$$         call wkonfchk(alat,plat,nbas,lmxax,lmxa,nsp,konft)
c$$$!!
c$$$         open(newunit=ifi,file='HAMindex0',form='unformatted')
c$$$         write(ifi) alat,plat,nbas,lmxax,nsp,ngrp,ndima,norb
c$$$         write(ifi) konft(0:lmxax,1:nbas,1:nsp),lmxa(1:nbas)
c$$$         write(ifi) (ssite(ib)%class,ib=1,nbas)
c$$$         write(ifi) iqnum(:,1:norb),strn4c(1:norb)
c$$$         write(ifi) (trim(sspec(ssite(ib)%spec)%name),ib=1,nbas)
c$$$         write(ifi) symops(1:3,1:3,1:ngrp),invgx(1:ngrp),shtvg(1:3,1:ngrp)
c$$$         close(ifi)
c$$$c$$$         open(newunit=ifisym,file='SYMOPS')
c$$$c$$$         write(ifisym,*) ngrp
c$$$c$$$         do ig = 1,ngrp
c$$$c$$$           write(ifisym,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
c$$$c$$$           do i=1,3
c$$$c$$$             write(ifisym,"(3e24.16)") symops(i,1:3,ig)
c$$$c$$$           enddo
c$$$c$$$         enddo
c$$$c$$$         close(ifisym)
c$$$         return
c$$$      endif
      
!---------------------------------------------------------
!!--- MTO part ---- obtain norbmto, lxx,kxx ----
      allocate( ltabx(n0*nkap0,nbas),ktabx(n0*nkap0,nbas),offlx(n0*nkap0,nbas))
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
      norbmto=0
      ndimham = 0 !dimension of mto part of hamiltonian
      allocate(offH(nbas+1)) !offH looks 
      offH=0
      do  ib = 1, nbas
        is=ssite(ib)%spec
c        spid(ib)=sspec(is)%name
        call orbl(ib,0,ldim,iprmb,norb,ltabx(:,ib),ktabx(:,ib),off,offlx(:,ib),ndim) !iprmb
        do  iorb = 1, norb !(ib,irob) specify a block of MTO part Hamiltonian
          norbmto=norbmto+1
          ibastab(norbmto)= ib
          ltab(norbmto)   = ltabx(iorb,ib) !angular momentum l of (ib,iorb) block
          ktab(norbmto)   = ktabx(iorb,ib) !radial index of (ib,iorb) block
          offl(norbmto)   = offlx(iorb,ib) !offset to (ib,iorb) block
          nini = ndimham+ 1 
          ndimham = ndimham+ 2*ltab(norbmto)+1
          ibasindex(nini:ndimham) = ib
c ib,ltab(norbmto),ktab(norbmto), offl(norbmto)+1,ndimham,trim(spid(ib))
        enddo
        offH(ib+1) = ndimham !'starting index'-1 of (ib) block
      enddo
      offH(nbas+1) = ndimham
c ... reverse maping of offset-index for hamiltonian
      allocate(offlrev(nbas,0:lxx,kxx))
      do iorb=1,norbmto
        ibas = ibastab(iorb)
        l   = ltab(iorb)
        k   = ktab(iorb)
        offlrev(ibas,l,k)= offl(iorb)
      enddo
!! Get rotation matrix dlmm.  We assume nl=lmxa+1.
      lxxa=nl-1
      allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      call rotdlmm(symops,ngrp, nl, dlmm)

!! Not GW mode ------------------------------------
      if(jobgw<0) return !for sigm mode, dlmm needed.
!---------------------------------------------------------
      
!! --- PW part. info for eigenfunctions are expanded as MTpart+PWpart.!feb2012takao
      inquire(file='QGpsi',EXIST=qpgexist)  !feb2012takao
      if(.not.qpgexist) then
        call writehamindex()
        goto 2001 
      endif
!! q on mesh and shortened q.
      open(newunit=ifiqg,file='QGpsi',form='unformatted')
      read(ifiqg) nqnum, ngpmx ,QpGcut_psi, nqbz___, nqi !,imx !,nqibz
!! we have two set of data for original qxx in QGpsi and their shortened.
      if(allocated(qq)) deallocate(qq)
      nqtt=nqnum !nqnum*2 !doubled. second series nqnum+1:2*nqnum are for shortened q.
      nkt=nqtt
      allocate( qtti(3,nqi), qq(3,nqtt),iqtt(nqtt) )
      iqi=0
      do  iq = 1, nqnum
        read(ifiqg)  qxx,ngp,irr  ! q, and number of G vectors for
        if(irr/=0) then
          iqi=iqi+1
          qtti(:,iqi)=qxx
          iqtt(iqi)=iq
        endif
        read(ifiqg)
        qq(:,iq)=qxx
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
        do iqi=1,nqi
          q=qtti(:,iqi)
          iqq=iqtt(iqi)
          iqi_=iqi
          do ig=1,ngrp
            call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
            if(sum(abs(qx))<tolq()) then
              igg=ig
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
      pwmode=ham_pwmode
      if(pwmode==0.or.pwemax<1d-8) then
        if(allocated(napwk)) deallocate(napwk)
        allocate(napwk(nkt))
        napwk=0
        napwmx=0
        if(master_mpi) then
           print *,'pwmode=0 writehamindex'
           call writehamindex() !sep2012takao
        endif   
        return
      endif
!! for APW rotation.  ! ... Get igv2(3,iapw,ikt). pwmode>=10 only
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
        if(pwmode<10) then
           ppin=matmul(transpose(plat),qq(:,ikt))
           call shortn3(ppin,noutmx, nout,nlatout)
        endif
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
        do ig  = 1,napwk( ikt )
          nnn  = igv2(1:3, ig, ikt)
          igv2rev( nnn(1), nnn(2),nnn(3), ikt) = ig
        enddo
      enddo
      if(master_mpi) call writehamindex()
      
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
        deallocate(symops, invgx,miat,tiat,shtvg,  dlmm )
!!---- get space group information ---- translation informations also in miat tiat invgx, shtvg
        ngrp_original=ngrp
        ngrp      = ngrpaf ! Overwrite ngrp by ngrpaf (>ngrp because we treat AF pairs are in
                           ! the same class.
        allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
        call mptauof ( symtmp , ngrp, plat , nbas , rv_a_opos , iclasstaf
     .   , miat , tiat , invgx , shtvg )
        if(master_mpi) then
          write(6,*)
          write(6,"(' ngrp for SYMGRP+GYMGRPAF, ngrp for SYMGRP=',2i5)") ngrp, ngrp_original
          do ig=1,ngrp
          write(6,"(a,i3,a,100i3)")'ig=',ig, ' miat=',miat(:,ig)
          enddo
        endif  
        lxxa=nl-1
        allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
        call rotdlmm(symtmp,ngrp, nl, dlmm) !rotation matrix dlmm.  We assume nl=lmxa+1.
      endif
      end subroutine m_hamindex_init
      
      
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
!> get index ikt such that for qin(:)=qq(:,ikt)
      integer function getikt(qin) !return
      intent(in)::            qin
      integer::i
      real(8):: qin(3)
      getikt=-99999
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

!> write info for wave rotation. (internal subroutine)
      subroutine writehamindex()
      integer(4):: ifi
      logical::pmton
      logical,save:: done=.false.
      if(done) call rx('writehamindex is already done')
      done=.true.
      zbak=zbak_read
      ndham=ndham_read
      open(newunit=ifi,file='HAMindex',form='unformatted')
      write(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto,pwmode,zbak,ndham
      write(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
      write(ifi)lxxa
      write(ifi)dlmm
      write(ifi)ibastab,ltab,ktab,offl,offlrev !for rotation of MTO. recovered sep2012 for EIBZ for hsfp0
      write(ifi)qq !,ngvecp,ngvecprev
      write(ifi)plat,qlat,napwmx
      if(napwmx/=0) then !for APW rotation used in rotwvigg
         write(ifi) igv2,napwk,igv2rev
      endif
!xxx      
      write(ifi) alat,rv_a_opos
      
      close(ifi)
      end subroutine writehamindex

!> Read info of PMT Hamiltoninan
      subroutine readhamindex()
!! == read info for wave rotation. feb2012takao ==
      integer(4):: ifi,nkt
      logical::pmton
      logical,save:: done=.false.
      if(done) call rx('readhamindex is already done')
      done=.true.
      open(newunit=ifi,file='HAMindex',form='unformatted')
      read(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto,pwmode,zbak,ndham
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
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine pvsug1(nbas,lmxa,ipc,ipcx)
C- Remakes class table, expunging classes with floating orbitals.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   lmxa  :augmentation l-cutoff
Ci   ipc   :class table: site ib belongs to class ipc(ib)
Co Inputs/Outputs
Co   ipcx  :expunged class table: classes with lmxa=-1 are expunged
Co         :and the remaining classes are sequentially renumbered
Co         :preserving the order of the remaining classes.
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Mar 07
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,lmxa(nbas),ipc(nbas),ipcx(nbas)
C ... Local parameters
      integer i,ip,ic,ipskip,nelim
      integer prm(nbas)
      call ivshel(1,nbas,ipc,prm,.true.)
      nelim = 0
      ipskip = 0
      do  i = 1, nbas
        ip = prm(i)+1
        ic = ipc(ip)
        if (lmxa(ip) .lt. 0) then
          if (ipskip .ne. ic) nelim = nelim+1
          ipskip = ic
          ipcx(ip) = -1
        else
          ipcx(ip) = ic - nelim
        endif
      enddo
      end subroutine pvsug1
!!
      
      end module m_hamindex


