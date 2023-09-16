module m_mksym !crystal symmetry data given by call m_mksym_init.  nbas (atomic sites)-> nspec (species) -> nclass (class)
  public :: m_mksym_init
  !NOTE: We have ngrp+ngrpAF symmetry for symops... ngrp for lattice, ngrpAF is extra symmetry for AF (spin-flip symmetry)
  integer,allocatable,protected   :: oics(:)    ! ispec= ics(iclass) gives spec for iclass.
  real(8),allocatable,protected   :: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:),dlmm(:,:,:,:)
  integer,allocatable,protected :: invgx(:),miat(:,:), oistab (:,:)   ! j= istab(i,ig): site i is mapped to site j by grp op ig
  integer,allocatable,protected   :: iclasst(:),ltab(:),ktab(:),offl(:), offlrev(:,:,:),ibastab(:)
  integer,protected::                nclasst !number of equivalent class
  integer,protected :: ngrp  !# of lattice symmetry 
  integer,protected :: npgrp !# of lattice symmetry with artificial inversion (as time-reversal symmetry) if lc=1 =>mkqp
  
  logical,protected:: AFmode
  integer,protected:: ngrpAF !SpinFlip symmetry: y=matmul(symopsF(:,:,ig),x)+ ag_af(:,ig), ig=ngrp+1,ngrp+ngrpAF.

  integer,protected :: lxx,kxx,norbmto,lxxa !oribtal index
  integer,allocatable,protected:: ib_table(:),k_table(:),l_table(:)
contains
  subroutine m_mksym_init(prgnam) 
    use m_mpitk,only:   master_mpi
    use m_lgunit,only:  stdo
    use m_lmfinit,only: nspec,nbas, sstrnsymg,symgaf,ips=>iv_a_oips,slabl,iantiferro,addinv,lmxax
    use m_lattic,only:  plat=>lat_plat,rv_a_opos
    use m_mksym_util,only: mksym
    use m_ftox
    implicit none
    integer,parameter::  ngmx = 48
    character,intent(in)::  prgnam*(*)
    integer:: ibas,lc,j,iprint,nclass,ngrpAll,k,npgrpAll
    integer,parameter::recln=511
    logical ::cmdopt0,ipr10=.false.
    character strn*(recln),strn2*(recln),outs(recln)
    real(8):: osymgr (3,3,ngmx),oag(3,ngmx)
    real(8),parameter:: tol=1d-8
    integer,allocatable:: iclasstAll(:)
    call tcn('m_mksym_init')
    if(.NOT.prgnam=='LMFA') ipr10= iprint()>10 !this is only for master
    strn = 'find'
    if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
    if(cmdopt0('--nosym') .OR. cmdopt0('--pdos') ) strn = ' '
    lc=merge(1,0,addinv) ! Add inversion to get sampling k points. phi*. When we have TR with keeping spin \sigma, psi_-k\sigm(r) = (psi_k\sigma(r))^* 
    lxxa=lmxax
    if(master_mpi) call pshpr(60)
    allocate(iclasst(nbas),oics(nbas),oistab(nbas,ngmx))
    if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice: ========start========================== '
    call mksym(lc,slabl,strn,ips, iclasst,nclasst,npgrp,ngrp,oag,osymgr,oics,oistab)
      if(ipr10) write(stdo,"(a)") 'SpaceGroupSym of Lattice: ========end =========================== '
    allocate(symops,source=osymgr)
    allocate(ag,source=oag)
    
    AFmodeBlock: block
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
      real(8):: osymgrAll(3,3,ngmx),oagAll(3,ngmx)
      integer:: oicsAll(nbas),oistabAll(nbas,ngmx),ipsAF(nbas),iga,igall,nclassAll,ig
      AFmode=len_trim(symgaf)>0 
      if(AFmode) then
         strn2=trim(strn)//' '//trim(symgaf)
         if(ipr10) then
            write(stdo,*)
            write(stdo,"(a)") 'Add SpaceGroupSym ops by AF symmetry===start========= '
            write(stdo,"(a)") 'AF: Antiferro mode: SPGGRAF  = '//trim(symgaf)
            write(stdo,"(a)") 'AF:                 SPGGR all= '//trim(strn2)
            write(stdo,"(a,2i3)")  ('AF:  ibas,AF=',j,iantiferro(j),j=1,nbas)
         endif
         ipsAF = ips
         do j=1,nbas
            do k=j,nbas
               if( iantiferro(j)+iantiferro(k)==0) then
                  ipsAF(k) = ipsAF(j) !to drive mksym for AF mode (Assuming AF pairs with the same spec).
                  exit
               endif
            enddo
         enddo
         if(master_mpi) call pshpr(60)
         allocate(iclasstAll(nbas))
         if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice+AF: ========start========================== '
         call mksym(lc,slabl,strn2,ipsAF, iclasstAll,nclassAll, npgrpAll,ngrpAll,oagAll,osymgrAll,oicsAll,oistabAll) !Big symmetry for lattice+AF
         ngrpAF=ngrpAll-ngrp
         if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice+AF: ========end========================== '
         iga=ngrp
         do igall=1,ngrpAll !Pick up symmetry by AF
            if(any([( sum(abs(symops(:,:,ig)-osymgrALL(:,:,igall)))+sum(abs(ag(:,iga)-oagAll(:,igall)))<tol,ig=1,ngrp)]) ) cycle
            iga=iga+1
            symops(:,:,iga)  = osymgrALL(:,:,igall)
            ag(:,iga)      = oagALL(:,igall)
            oistab(:,iga)= oistabAll(:,igall)
         enddo
         if(iga/=ngrpAF+ngrp) call rx('ngrpAF+nggp/=ngrpAll')
      else
         ngrpAll=ngrp
         ngrpAF=0
      endif
    endblock AFmodeBlock
    MiatTiatDlmm:block
      allocate(miat(nbas,ngrpAll),tiat(3,nbas,ngrpAll),invgx(ngrpAll),shtvg(3,ngrpAll),dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa,ngrpAll))
      call            mptauof(symops,ngrp,  plat,nbas,rv_a_opos,iclasst,   miat,           tiat,             invgx,         shtvg) !Get spg info
      if(AFmode) call mptauof(symops,ngrpAF,plat,nbas,rv_a_opos,iclasstAll,miat(:,ngrp+1:),tiat(:,:,ngrp+1:),invgx(ngrp+1:),&
           shtvg(:,ngrp+1)) !Get spg info
!      write(6,*)'mmmmmmmmmmmmmmkkkkk111',invgx(1:ngrpAll)
!      write(6,*)'mmmmmmmmmmmmmmkkkkk222',ngrpAF,ngrpAll
      call rotdlmm(symops,ngrpAll, lxxa+1, dlmm) !!! Get rotation matrix Dlmm in real spherical harmonics.  !for sigm mode, dlmm needed.
      if(ipr10) write(stdo,*)
    endblock MiatTiatDlmm
    Orbitalindex: block! A block contains 2*l+1 orbitals. A block specified by (ibas,k,l) !k=1,2,3 is for EH,EH2,PZ
      use m_lmfinit,only: ldim=>nlmto, nbas,ispec,norbx,ltabx,ktabx,offlx,nl
      integer:: ib,iorb,is,k,l
      norbmto=0
      kxx=-1
      lxx=-1 
      do  ib = 1, nbas
         is=ispec(ib) 
         do iorb = 1, norbx(ib)
            norbmto = norbmto+1
            if(ltabx(iorb,ib)>lxx)  lxx = ltabx(iorb,ib)
            if(ktabx(iorb,ib)>kxx)  kxx = ktabx(iorb,ib) 
         enddo
      enddo
      !!--- make index table :norbmto is the total number of different type of MTOs !      allocate( ibasindex(ndimham))
      allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
      norbmto=0 !      ndimham = 0 !dimension of mto part of hamiltonian!      allocate(offH(nbas+1)) !offH looks!      offH=0
      do  ib = 1, nbas
         is=ispec(ib) 
         do  iorb = 1, norbx(ib) !(ib,irob) specify a block of MTO part Hamiltonian
            norbmto=norbmto+1
            ibastab(norbmto)= ib
            ltab(norbmto)   = ltabx(iorb,ib) !angular momentum l of (ib,iorb) block
            ktab(norbmto)   = ktabx(iorb,ib) !radial index of (ib,iorb) block
            offl(norbmto)   = offlx(iorb,ib) !offset to (ib,iorb) block
         enddo
      enddo
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
         ib_table(offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ib
         l_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ltab(iorb)
         k_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ktab(iorb)
      enddo
    endblock Orbitalindex
    call tcx('m_mksym_init')
  end subroutine m_mksym_init
end module m_mksym

