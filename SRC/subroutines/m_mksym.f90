module m_mksym !crystal symmetry data given by call m_mksym_init.  nbas (atomic sites)-> nspec (species) -> nclass (class)
  public :: m_mksym_init
  integer,allocatable,protected   :: iv_a_oics(:)        ! ispec= ics(iclass) gives spec for iclass.
  real(8),allocatable,protected   :: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:),dlmm(:,:,:,:),rv_a_osymgr (:,:,:)
  real(8) , allocatable,protected :: rv_a_oag (:,:)      ! translation for 1,...nsgrp
  integer , allocatable,protected :: invgx(:),miat(:,:),iv_a_oistab (:,:)   ! j= istab(i,ig): site i is mapped to site j by grp op ig
  integer,allocatable,protected   :: iclasst(:),ltab(:),ktab(:),offl(:), offlrev(:,:,:),ibastab(:)
  integer,protected:: lat_nsgrp,               ctrl_nclass,                  lat_npgrp
  !                   # of lattice symmetry   !number of equivalence class  !number of point group operation with adding iversion. Used in mkqp.f90
  integer,protected :: ngrp
  
  !antiferro sym y = matmul(sympsAF(:,:,ig),x)+ ag_af(:,ig), ig=1,ngrpAF
  logical:: AFmode
  integer,allocatable,protected:: iclasstAF(:) !=== AntiFerro class information 
  real(8),allocatable,protected:: symopsAF(:,:,:), agAF(:,:) 
  integer,protected:: ngrpAF       !antiferro symopsAF (For NiO, ngrpAF=24=2*lat_nsgrp
  real(8),allocatable:: dlmmAF(:,:,:,:)

  integer,protected :: lxx,kxx,norbmto,lxxa
  integer,allocatable:: ib_table(:),k_table(:),l_table(:)
contains
  subroutine m_mksym_init(prgnam) 
    use m_mpitk,only:   master_mpi
    use m_lgunit,only:  stdo
    use m_lmfinit,only: nspec,nbas, sstrnsymg,symgaf,iv_a_oips,slabl,iantiferro,addinv !input
    use m_lattic,only:  rv_a_opos !input
    use m_mksym_util,only: mksym
    use m_ftox
    ! note: Antiferro symmetric self-energy mode.
    !! To obtain self-consistency, it may be useful to keep AF condition during iteration.
    !! We need to set SYMGRPAF (still developing...)
    ! sample for nio
    ! SYMGRPAF i:(1,1,1) ! translation + inversion makes antiferro symmetry.
    ! SYMGRP r3d   !this keeps spin axis. (probably larger symmetry allowed for SO=0)
    ! STRUC  ALAT={a} PLAT= 0.5 0.5 1.0  0.5 1.0 0.5  1.0 0.5 0.5
    !        NBAS= 4  NSPEC=3
    ! SITE   ATOM=Niup POS=  .0   .0   .0   AF=1    <--- AF symmetric pair
    !        ATOM=Nidn POS= 1.0  1.0  1.0   AF=-1   <--- 
    !        ATOM=Mnup POS=  .0   .0   .0   AF=2   (if another AF pairs)
    !        ATOM=Mndn POS= 1.0  1.0  1.0   AF=-2 
    !        ATOM=O POS=  .5   .5   .5
    !        ATOM=O POS= 1.5  1.5  1.5
    !------------------------------------------------------------
    implicit none
    integer,parameter::  ngmx = 48
    character,intent(in)::  prgnam*(*)
    integer:: ibas,lc,j,iprint,nclass,npgrp,nsgrp
    integer,parameter::recln=511
    logical ::cmdopt0,ipr10=.false.
    character strn*(recln),strn2*(recln),outs(recln)
    call tcn('m_mksym_init')
    strn = 'find'
    if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
    if (cmdopt0('--nosym') .OR. cmdopt0('--pdos') ) strn = ' '
    lc=merge(1,0,addinv) ! Add inversion to get sampling k points. phi*. When we have TR with keeping spin \sigma, psi_-k\sigm(r) = (psi_k\sigma(r))^* 
    if(.NOT.prgnam=='LMFA') ipr10= iprint()>10 !this is only for master
    AFmode=len_trim(symgaf)>0 
    allocate( rv_a_oag(3,ngmx)      )
    allocate( rv_a_osymgr(3,3,ngmx) )
    allocate( iv_a_oics(nbas)  )
    allocate( iv_a_oistab(nbas,ngmx))
    AFsymblock :block
      integer:: inumaf,iv_a_oipsAF(nbas),j,k, ib
      if(AFmode) then
         strn2=trim(strn)//' '//trim(symgaf)
         iv_a_oipsAF=iv_a_oips
         if(ipr10) then
            write(stdo,*)
            write(stdo,"(a)") 'SpaceGroupSym including AF symmetry===start========= '
            write(stdo,"(a)") 'AF: ======== AF symmetry ===================== '
            write(stdo,"(a)") 'AF: Antiferro mode: SPGGRAF  = '//trim(symgaf)
            write(stdo,"(a)") 'AF:                 SPGGR all= '//trim(strn2)
            write(stdo,"(a,2i3)")  ('AF:  ibas,AF=',j,iantiferro(j),j=1,nbas)
         endif
         !NOTE: a little confusing. ! Module variables are written by mksymaf but overwitten by next call of mksym.
         inumaf = 0
         do j=1,nbas
            do k=j,nbas
               if( iantiferro(j)+iantiferro(k)==0) then
                  iv_a_oipsAF(k) = iv_a_oipsAF(j) !to drive mksym for AF mode (Assuming AF pairs with the same spec).
                  inumaf=inumaf+1
                  exit
               endif
            enddo
         enddo
         if(master_mpi) call pshpr(60)
         allocate(iclasstAF(nbas))
         call mksym(lc,slabl,strn2,iv_a_oipsAF, iclasstAF,nclass,npgrp,ngrpAF,rv_a_oag,rv_a_osymgr,iv_a_oics,iv_a_oistab)!nclass... are overwritten by next mksym
         allocate(symopsAF(3,3,ngrpAF),agAF(3,ngrpAF)) !ngrpAF  = number of lattice symmetry  + af symmetry.
         if(master_mpi) call poppr()
         if(master_mpi) write(stdo,"(a)")'AF: mksym, generator= SYMGRP+SYMGRPAF= '//trim(strn2)
         call dcopy ( ngrpAF * 9 , rv_a_osymgr , 1 , symopsAF , 1 )
         call dcopy ( ngrpAF * 3 , rv_a_oag ,    1 , agAF , 1 )
         if(master_mpi) write(stdo,"(a,i3)")'AF: ngrpaf=',ngrpAF
         if(ipr10) write(stdo,"(a)") 'SpaceGroupSym of AF:      ========end =========================== '
         if(ipr10) write(stdo,"(a)")
      endif
    endblock AFsymblock
    if(master_mpi) call pshpr(60)
    allocate(iclasst(nbas))
    if(ipr10) write(stdo,"(a)")  'SpaceGroupSym of Lattice: ========start========================== '
    call mksym(lc,slabl,strn,iv_a_oips, iclasst,ctrl_nclass,lat_npgrp,lat_nsgrp, rv_a_oag,rv_a_osymgr,iv_a_oics,iv_a_oistab)
!    ctrl_nclass=nclass
    if(ipr10) write(stdo,"(a)") 'SpaceGroupSym of Lattice: ========end =========================== '
    if(ipr10) write(stdo,*)
    if(master_mpi) call poppr
    !------------------------
    mksymblock: block
!      use NaNum,only: NaN       !for initialization, but not working well
      use m_lmfinit,only: ldim=>nlmto, nbas,ispec,norbx,ltabx,ktabx,offlx,nl
      use m_lattic,only: plat=>lat_plat,rv_a_opos
      integer:: ib,iorb,is,k,l!,ndimham
!      use m_suham,only: ndham_read=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
!      integer:: ngrp_original!,pwmode,ndham !ngrpaf,
!      integer,allocatable::  offH (:)
!      real(8),allocatable:: symtmp(:,:,:)
!      real(8),allocatable:: symops_af(:,:,:), ag_af(:,:)
!      integer:: nqi=NaN, nqnum=NaN, nqtt=NaN, ndimham=NaN, napwmx=NaN,ngpmx=NaN
!           iqimap(:),iqmap(:),igmap(:),ibasindex(:), &
!           igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)! for rotation of evec       
!           qq(:,:), qtt(:,:),qtti(:,:) 
!      logical:: readhamindex_init=.false.
!      logical:: debug=.false.
!      integer:: ibas,k,l,ndim,ipr,nglob,off,offs,specw,fieldw,iorb,offsi,ib,is, norb
!      integer:: nkabc(3),nkp,lshft(3),napwx,ig,nini,nk1,nk2,nk3,ik1,ik2,ik3,ikt,nkt
!      integer:: ifi,ifisym,i,ifiqibz,igg,iqq,iqi,irr,iqi_ !,jobgw
!      integer:: iout,iapw ,iprint,ngadd,igadd!,igaf !,nout,nlatout(3,noutmx)
!      integer:: ngp, ifiqg,iq,nnn(3),ixx,ndummy,nqbz___ ,ifatomlist
!      integer,allocatable:: iqtt(:), kv(:)
!      real(8):: pwgmax, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)!pwgmin=0d0, 
!      real(8):: dum,qb(3,3),ddd(3),ppin(3), rlatp(3,3),xmx2(3),qqq(3),diffs,ddf
!      logical:: siginit, qpgexist, llmfgw,prpushed
!      character(8)::  spid(nbas)
!      integer::ndima,lmxax,npqn,ificlass,nat,lmaxa,ipqn,ifinlaindx,isp,konf
!    plat=lat_plat
!    qlat=lat_qlat
!    nsp=nsp_in
      ngrp=lat_nsgrp !note nsgrp given in mksym.F is without inversion.
      allocate(symops(3,3,ngrp),ag(3,ngrp))
      call dcopy ( ngrp * 9 , rv_a_osymgr , 1 , symops , 1 )
      call dcopy ( ngrp * 3 , rv_a_oag , 1 , ag , 1 )
      allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp)) !iclasst(nbas),
      !    do ib=1,nbas !       iclasst(ib)=ssite(ib)%class!    enddo
      !! get space group information ---- translation informations also in miat tiat invgx, shtvg
      call mptauof(symops,ngrp,plat,nbas, rv_a_opos, iclasst, miat , tiat , invgx , shtvg )
      norbmto=0
      kxx=-1
      lxx=-1
!      ndimham = 0               !dimension of mto part of hamiltonian
      do  ib = 1, nbas
         is=ispec(ib) 
         do iorb = 1, norbx(ib)
            norbmto = norbmto+1
            if(ltabx(iorb,ib)>lxx)  lxx = ltabx(iorb,ib)
            if(ktabx(iorb,ib)>kxx)  kxx = ktabx(iorb,ib)
!            ndimham = ndimham+ 2*ltabx(iorb,ib)+1
         enddo
      enddo
      !!--- make index table :norbmto is the total number of different type of MTOs
!      allocate( ibasindex(ndimham))
      allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
      norbmto=0
!      ndimham = 0 !dimension of mto part of hamiltonian
!      allocate(offH(nbas+1)) !offH looks
!      offH=0
      do  ib = 1, nbas
         is=ispec(ib) 
         do  iorb = 1, norbx(ib) !(ib,irob) specify a block of MTO part Hamiltonian
            norbmto=norbmto+1
            ibastab(norbmto)= ib
            ltab(norbmto)   = ltabx(iorb,ib) !angular momentum l of (ib,iorb) block
            ktab(norbmto)   = ktabx(iorb,ib) !radial index of (ib,iorb) block
            offl(norbmto)   = offlx(iorb,ib) !offset to (ib,iorb) block
!            nini = ndimham+ 1
!            ndimham = ndimham+ 2*ltab(norbmto)+1
!            ibasindex(nini:ndimham) = ib           ! ib,ltab(norbmto),ktab(norbmto), offl(norbmto)+1,ndimham,trim(spid(ib))
         enddo
!         offH(ib+1) = ndimham !'starting index'-1 of (ib) block
      enddo
!      offH(nbas+1) = ndimham
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
!         spid(ib)=slabl(is) 
         ib_table(offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ib
         l_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ltab(iorb)
         k_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ktab(iorb)
      enddo
      lxxa=nl-1
      allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      call rotdlmm(symops,ngrp, nl, dlmm) !!! Get rotation matrix dlmm.  We assume nl=lmxa+1. !for sigm mode, dlmm needed.
      !
      AFsymPart: if(AFmode) then !Symmetry for AF for GW. Order AF symmetry operation after normal one. jun2015takao
!         igadd=ngrp
!         do igaf=1,ngrpaf
!            do ig=1,ngrp
!               diffs=sum(abs(symops_af(:,:,igaf)-symops(:,:,ig)))
!               if(diffs<1d-6) goto 1013
!            enddo
!            igadd=igadd+1
!            symtmp(:,:,igadd)=symops_af(:,:,igaf)
!1013        continue
!         enddo
!         if(igadd/=ngrpaf) call rx('m_hamindex: strange. bug igadd/=ngrpaf')
!         if(master_mpi) write(stdo,*) '-----SYMGRPAF mode ---- # of additional symmetry=',igadd
         allocate( dlmmAF( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrpAF))
         call rotdlmm(symopsAF,ngrpAF, nl, dlmmAF) !rotation matrix dlmmAF.  We assume nl=lmxa+1.
      endif AFsymPart
      
    endblock mksymblock
    call tcx('m_mksym_init')
  end subroutine m_mksym_init
end module m_mksym

