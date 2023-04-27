module m_hamindex   !hamiltonian index read/write for successive GW calculaitons.
  use m_lmfinit,only: ham_pwmode,pwemax,ldim=>nlmto,noutmx,nsp_in=>nsp,stdo, &
       alat=>lat_alat,nl,ctrl_nbas=>nbas,ispec,sspec=>v_sspec,n0,nkap0,zbak_read=>zbak,slabl,z
  use m_lattic,only: lat_qlat,lat_plat,rv_a_opos
  use NaNum,only: NaN       !for initialization, but not working well
  use m_suham,only: ndham_read=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
  use m_lmfinit,only:norbx,ltabx,ktabx,offlx
  use m_lgunit,only:stdo
  use m_ftox
  public:: m_hamindex_init, Readhamindex, getikt
  integer,protected,allocatable,public:: ib_table(:),k_table(:),l_table(:)
  integer,protected,public:: ngrpaf,ngrp_original,pwmode,ndham
  integer,protected,public:: nqi=NaN, nqnum=NaN, ngrp=NaN, lxx=NaN, kxx=NaN,norbmto=NaN, &
       nqtt=NaN, ndimham=NaN, napwmx=NaN, lxxa=NaN, ngpmx=NaN, imx=NaN,nbas=NaN
  integer,allocatable,protected,public:: iclasstaf(:), offH (:), &
       ltab(:),ktab(:),offl(:), offlrev(:,:,:),ibastab(:), & !iclasst(:),
       iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:), &
       igv2(:,:,:),napwk(:),igv2rev(:,:,:,:),igvx(:,:)
  real(8),allocatable,protected,public:: symops_af(:,:,:), ag_af(:,:), &
       symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:), qtt(:,:),qtti(:,:)
  real(8),protected,public:: plat(3,3)=NaN,qlat(3,3)=NaN,zbak
  logical,protected,public:: readhamindex_init=.false.
  private
  logical,private:: debug=.false.
contains
  subroutine m_hamindex_init(jobgw)
    use m_mksym,only: rv_a_osymgr,rv_a_oag,lat_nsgrp, iclasstaf_,symops_af_,ag_af_,ngrpaf_,iclasst
    use m_struc_def
    use m_MPItk,only: master_mpi
    !!-- Set up m_hamiltonian. Index for Hamiltonian. --
    !!  Generated index are stored into m_hamindex 
    !!  Only include q-point information for GW (QGpsi).
    !!#Inputs
    !r As you see in subroutine rotwvigg, the index for Hamiltonian reads as;
    !      do iorb=1,norbmto             !orbital-blocks are specified by ibas, l, and k.
    !        ibas  = ibastab(iorb)
    !         l    = ltab(iorb)
    !         k    = ktab(iorb)        !kappa (diffent MTO index for each l)
    !        init1 = offl(iorb)+1      !starting index for the block iorb
    !        iend1 = offl(iorb)+2*l+1  !end of the block for the iorb
    !      enddo
    implicit none
    integer:: ibas,k,l,ndim,ipr,nglob,off,offs,specw,fieldw,iorb,offsi,ib,is, norb,nsp
    integer:: nkabc(3),nkp,lshft(3),napwx,ig,nini,nk1,nk2,nk3,ik1,ik2,ik3,ikt,nkt
    integer:: i_copy_size,i_spacks,i_spackv,ifi,nbas_in,ifisym,i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw
    integer:: iout,iapw ,iprint,ngadd,igadd,igaf !,nout,nlatout(3,noutmx)
    integer:: ngp, ifiqg,iq,nnn(3),ixx,ndummy,nqbz___ ,ifatomlist
    integer,allocatable:: iqtt(:), kv(:)!ltabx(:,:),ktabx(:,:),offlx(:,:),
    real(8):: pwgmax, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)!pwgmin=0d0, 
    real(8):: dum,qb(3,3),ddd(3),ppin(3), rlatp(3,3),xmx2(3),qqq(3),diffs,ddf
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
    call tcn('m_hamindex_init')
    nbas=ctrl_nbas !nbas_in
    ngrp=lat_nsgrp !note nsgrp given in mksym.F is without inversion.
    plat=lat_plat
    qlat=lat_qlat
    nsp=nsp_in
    allocate(symops(3,3,ngrp),ag(3,ngrp))
    call dcopy ( ngrp * 9 , rv_a_osymgr , 1 , symops , 1 )
    call dcopy ( ngrp * 3 , rv_a_oag , 1 , ag , 1 )
    allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp)) !iclasst(nbas),
!    do ib=1,nbas !       iclasst(ib)=ssite(ib)%class!    enddo
    !! get space group information ---- translation informations also in miat tiat invgx, shtvg
    call mptauof(symops, ngrp , plat , nbas , rv_a_opos, iclasst, miat , tiat , invgx , shtvg )
    norbmto=0
    kxx=-1
    lxx=-1
    ndimham = 0               !dimension of mto part of hamiltonian
    do  ib = 1, nbas
       is=ispec(ib) 
       do iorb = 1, norbx(ib)
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
       is=ispec(ib) 
       do  iorb = 1, norbx(ib) !(ib,irob) specify a block of MTO part Hamiltonian
          norbmto=norbmto+1
          ibastab(norbmto)= ib
          ltab(norbmto)   = ltabx(iorb,ib) !angular momentum l of (ib,iorb) block
          ktab(norbmto)   = ktabx(iorb,ib) !radial index of (ib,iorb) block
          offl(norbmto)   = offlx(iorb,ib) !offset to (ib,iorb) block
          nini = ndimham+ 1
          ndimham = ndimham+ 2*ltab(norbmto)+1
          ibasindex(nini:ndimham) = ib           ! ib,ltab(norbmto),ktab(norbmto), offl(norbmto)+1,ndimham,trim(spid(ib))
       enddo
       offH(ib+1) = ndimham !'starting index'-1 of (ib) block
    enddo
    offH(nbas+1) = ndimham
    allocate(offlrev(nbas,0:lxx,kxx))
    do iorb=1,norbmto ! ... reverse maping of offset-index for hamiltonian
       ibas = ibastab(iorb)
       l   = ltab(iorb)
       k   = ktab(iorb)
       offlrev(ibas,l,k)= offl(iorb)
    enddo
    allocate(ib_table(ldim),l_table(ldim),k_table(ldim))
    do iorb = 1, norbmto      !Total number of MTO's (without m)
       ib   = ibastab(iorb)
       is   = ispec(ib) 
       spid=slabl(is) 
       ib_table(offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ib
       l_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ltab(iorb)
       k_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ktab(iorb)
    enddo
    lxxa=nl-1
    allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
    call rotdlmm(symops,ngrp, nl, dlmm) !!! Get rotation matrix dlmm.  We assume nl=lmxa+1. !for sigm mode, dlmm needed.
    if(jobgw<0) then ! Not GW mode 
       call tcx('m_hamindex_init')
       return                 
    endif
    WriteHamindexBlock: block
      real(8):: tolq
      inquire(file='QGpsi',EXIST=qpgexist) ! ------------ GW mode   
      if( .NOT. qpgexist) goto 2001 ! skip writehamindex
      open(newunit=ifiqg,file='QGpsi',form='unformatted') ! q on mesh and shortened q.
      read(ifiqg) nqnum, ngpmx ,QpGcut_psi, nqbz___, nqi
      if(allocated(qq)) deallocate(qq)
      nqtt=nqnum 
      nkt=nqtt
      allocate( qtti(3,nqi), qq(3,nqtt),iqtt(nqtt) )
      iqi=0
      do  iq = 1, nqnum
         read(ifiqg)  qxx,ngp,irr  ! read q and number of G vectors (irr=1 meand irreducible points)
         if(irr/=0) then
            iqi=iqi+1
            qtti(:,iqi)=qxx
            iqtt(iqi)=iq
         endif
         read(ifiqg)
         qq(:,iq)=qxx
         if(master_mpi)write(stdo,"(' qq=',i5,3f10.5)") iq,qq(:,iq)
      enddo
      close(ifiqg)
      allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
      platt= transpose(plat) !inverse of qlat
      allocate(qtt(3,nqtt))
      qtt(:,1:nqtt)=qq(:,1:nqtt)
      do i=1,nqtt !Generate info for rotwv and write 
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
                  if(ddf-nint(ddf)>1d-8)write(stdo,ftox)'qxqx',ftof(q),ftof(qtarget),ftof(qtarget-matmul(symops(:,:,ig),q))
                  goto 2012
               endif
            enddo
         enddo
         if(master_mpi) then !error exit
            write(stdo,"(a,3f7.3)")'gen_ham: qtarget cannot found. Need to add SYMGRP explicitly (for SO=1), '//&
                 'or You have to delete inconsistent QGpsi. qtarget=',qtarget
            write(stdo,*)'gen_hamindex: qtarget can not found by SYMOPS.'
            write(stdo,"('qq20 ',3d16.8,2x,3d16.8,2x,3d16.8)")q,qtarget,matmul(platt,(qtarget-matmul(symops(:,:,ig),q)))
            do ig=1,ngrp
               call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
               write(stdo,"('qqqq2 ',3d16.8,2x,3d16.8,2x,3d16.8)") qtarget-matmul(symops(:,:,ig),q),qx
            enddo
            call rx('gen_hamindex: you may need to repeat echo 1|qg4gw, when you changed SYMOPS.')
         endif
2012     continue
         iqmap(i)=iqq
         iqimap(i)=iqi_
         igmap(i)=igg
      enddo

      !! === rotation of APW. (not the IPW part for GW).===
      pwmode=ham_pwmode
      if(master_mpi) call pshpr(0) !print index is pushed to be zero
      if(mod(pwmode,10)==0 .OR. pwemax<1d-8) then
         if(allocated(napwk)) deallocate(napwk)
         allocate(napwk(nkt))
         napwk=0
         napwmx=0
      else ! for APW rotation.  ! ... Get igv2(3,iapw,ikt). pwmode>=10 only
         if(master_mpi) print *,' gen_hamindex goto APW part: pwmode pwemax=',pwmode,pwemax !pwemin
         if(allocated(napwk)) deallocate(napwk,igv2,igv2rev)
         allocate( napwk(nkt))
         pwgmax = pwemax**.5
         do ikt=1,nkt
            qqq = merge(qq(:,ikt), 0d0, mod(pwmode/10,10)==1)
            call getgv2(alat,plat,qlat,qqq, pwgmax,1, napwk(ikt),dum)
            !call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,0, napwk(ikt),dum,dum,dum)!,dum)
         enddo
         napwmx=maxval(napwk)
         allocate( igv2(3,napwmx,nkt))!, kv(3*napwmx),igvx(napwmx,3))
         do ikt = 1,nkt
            qqq = merge(qq(:,ikt),0d0,mod(pwmode/10,10) == 1)
            !igvx=0
            call getgv2(alat,plat,qlat,qqq, pwgmax,2, napwk(ikt),igv2(:,:,ikt)) 
            !call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,1,napwmx,napwk(ikt),kv,dum,igvx) 
            !igv2(:,:,ikt) = transpose(igvx)
         enddo
!         deallocate(kv,igvx)
         ! ... 
         imx=-999
         do ikt = 1,nkt
            ixx = maxval( abs(igv2(1:3,1:napwk(ikt),ikt)))
            if(ixx>imx) imx=ixx
         enddo
         allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt),source=999999 ) !Reverse table of igv2 --->igv2rev
         do ikt = 1,nkt
            do ig  = 1,napwk( ikt )
               nnn  = igv2(1:3, ig, ikt)
               igv2rev( nnn(1), nnn(2),nnn(3), ikt) = ig
            enddo
         enddo
      endif
      if(master_mpi) call writehamindex()
2001  continue
      if(master_mpi) call poppr !print index is poped.
    endblock WriteHamindexBlock
    AFsymPart: if(allocated(symops_af)) then !Symmetry for AF for GW. Order AF symmetry operation after normal one. jun2015takao
       ! Caution: symops,and so on are overwritten. 
       !          writehamindex() already wrote HAMindex which is just for SYMGRP.
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
             if(diffs<1d-6) goto 1013
          enddo
          igadd=igadd+1
          symtmp(:,:,igadd)=symops_af(:,:,igaf)
1013      continue
       enddo
       if(igadd/=ngrpaf) call rx('suham: strange. bug igadd/=ngrpaf')
       if(master_mpi) write(stdo,*) '-----SYMGRPAF mode ---- # of additional symmetry=',igadd
       deallocate(symops, invgx,miat,tiat,shtvg,  dlmm )
       ngrp_original=ngrp ! get space group information ---- translation informations also in miat tiat invgx, shtvg
       ngrp      = ngrpaf ! Overwrite ngrp by ngrpaf (>ngrp because we treat AF pairs are in the same class.
       allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
       call mptauof ( symtmp , ngrp, plat , nbas , rv_a_opos , iclasstaf, miat , tiat , invgx , shtvg )
       if(master_mpi) then
          write(stdo,*)
          write(stdo,ftox)' ngrp for SYMGRP+GYMGRPAF=',ngrp,'ngrp for SYMGRP=',ngrp_original
          do ig=1,ngrp
             write(stdo,ftox)'ig=',ig, ' miat=',miat(:,ig)
          enddo
       endif
       lxxa=nl-1
       allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
       call rotdlmm(symtmp,ngrp, nl, dlmm) !rotation matrix dlmm.  We assume nl=lmxa+1.
    endif AFsymPart
    call tcx('m_hamindex_init')
  end subroutine m_hamindex_init
  integer function getikt(qin) !return !> get index ikt such that for qin(:)=qq(:,ikt)
    intent(in)::          qin
    integer::i
    real(8):: qin(3)
    getikt=-99999
    do i=1, nqnum !*2 !nkt
       if(debug) write(stdo,*)i,qin, qq(:,i)
       if(sum (abs(qin-qq(:,i)))<1d-8) then
          getikt=i
          return
       endif
    enddo
    write(stdo,*)' getikt: xxx error nqnum qin=',nqnum,qin
    do i=1, nqnum !*2 !nkt
       write(*,"('i qq=',i3,3f11.5)")i, qq(:,i)
    enddo
    call rx( ' getikt can not find ikt for given q')
  end function getikt
  subroutine writehamindex() !write info for wave rotation. (internal subroutine)
    integer:: ifi
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
    write(ifi) alat,rv_a_opos
    close(ifi)
  end subroutine writehamindex
  subroutine readhamindex() !Read info of PMT Hamiltoninan for wave rotation. 
    integer:: ifi,nkt
    logical::pmton
    logical,save:: done=.false.
    if(done) call rx('readhamindex is already done')
    done=.true.
    readhamindex_init=.true.
    open(newunit=ifi,file='HAMindex',form='unformatted')
    read(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto,pwmode,zbak,ndham
    allocate(symops(3,3,ngrp),ag(3,ngrp),qtt(3,nqtt),qtti(3,nqi))
    allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
    allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
    write(stdo,ftox) 'ngrp=',ngrp
    read(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
    allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
    allocate( offlrev(nbas,0:lxx,kxx))
    read(ifi) lxxa
    allocate( dlmm(-lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
    read(ifi) dlmm
    read(ifi)ibastab,ltab,ktab,offl,offlrev
    !      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
    !      allocate( ngvecp(3,ngpmx,nqnum) )
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
end module m_hamindex


