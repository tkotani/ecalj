!>Crystal symmetry data are stored by call m_mksym_init. NOTE:nbas (atomic sites)-> nspec (species) -> nclass (class)
module m_mksym 
  use m_cmdopt_registry, only: c0_nosym, c0_pdos
  public :: m_mksym_init
  private:: mksym
  integer,allocatable,protected :: oics(:)    ! ispec= ics(iclass) gives spec for iclass.
  real(8),allocatable,protected :: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:),dlmm(:,:,:,:)
  integer,allocatable,protected :: invgx(:),miat(:,:), oistab (:,:)   ! j= istab(i,ig): site i is mapped to site j by grp op ig
  integer,allocatable,protected :: iclasst(:)
  integer,protected::              nclasst !number of equivalent class
  integer,protected :: ngrp  !# of lattice symmetry 
  integer,protected :: npgrp !# of lattice symmetry with artificial inversion (as time-reversal symmetry) if lc=1 =>mkqp
  logical,protected:: AFmode
  integer,protected:: ngrpAF
  !NOTE: We have ngrp+ngrpAF symmetry for symops... ngrp for lattice, ngrpAF is extra symmetry for AF (spin-flip symmetry)
  !  ngrp+ngrpAF is the all symmetry for SYMGRP+SYMGRPAF  (SYMGRPAF needs AF index for AFpairs
  !  SpinFlip symmetry: y=matmul(symopsF(:,:,ig),x)+ ag_af(:,ig), ig=ngrp+1,ngrp+ngrpAF.
  !NOTE for Antiferro symmetric self-energy mode.
  !To obtain self-consistency, it may be useful to keep AF condition during iteration. We need to set SYMGRPAF.
  ! This is a sample for NiO
  ! SYMGRPAF i:(1,1,1) ! translation + inversion for spin-flip symmetry.
  ! SYMGRP r3d          ! this keeps spin axis. (probably larger symmetry allowed for SO=0)
  ! STRUC  ALAT={a} PLAT= 0.5 0.5 1.0  0.5 1.0 0.5  1.0 0.5 0.5
  !        NBAS= 4  NSPEC=3
  ! SITE   ATOM=Niup POS=  .0   .0   .0   AF=1    <--- AF symmetric pair
  !        ATOM=Nidn POS= 1.0  1.0  1.0   AF=-1   <--- 
  !        ATOM=Mnup POS=  .0   .0   .0   AF=2   (if another AF pairs)
  !        ATOM=Mndn POS= 1.0  1.0  1.0   AF=-2 
  !        ATOM=O POS=  .5   .5   .5
  !        ATOM=O POS= 1.5  1.5  1.5
contains
  subroutine m_mksym_init() 
    use m_mpi,only:   master_mpi
    use m_lgunit,only:  stdo
    use m_lmfinit,only: nspec,nbas, sstrnsymg,symgaf,ips=>iv_a_oips,slabl,iantiferro,addinv,lmxax
    use m_lattic,only:  plat=>lat_plat,rv_a_opos
    use m_symderive,only: mptauof,rotdlmm
    use m_ftox
    implicit none
    integer:: ngmx ! size of the operation arrays: 48 rotations x (pure translations <= nbas), x 2 for the added inversion (2026-10-02 05:46, step S4)
!    character,intent(in)::  prgnam*(*)
    integer:: ibas,lc,j,iprint,nclass,ngrpTotal,k,npgrpAll
    integer,parameter::recln=511
    logical ::ipr10=.false.
    character strn*(recln),strn2*(recln),outs(recln)
    real(8),allocatable:: osymgr(:,:,:), oag(:,:)
    real(8),parameter:: tol=1d-4
    integer,allocatable:: iclasstAll(:)
    call tcn('m_mksym_init')
    ipr10= iprint()>10 
    strn = 'find'
    if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
    if(c0_nosym .OR. c0_pdos ) strn = ' '
    lc=merge(1,0,addinv) ! Add inversion to get sampling k points. phi*. When we have TR with keeping spin \sigma, psi_-k\sigm(r) = (psi_k\sigma(r))^* 
    !lmxax=lmxax
    if(master_mpi) call pshpr(60)
    ngmx = 96*nbas
    allocate(osymgr(3,3,ngmx), oag(3,ngmx))
    allocate(iclasst(nbas),oics(nbas),oistab(nbas,ngmx))
    if(ipr10) write(stdo,"(a)")  'SpaceGroupSym of Lattice: ========start========================== '
    if(ipr10) write(stdo,"(a)") ' SYMGRP = '//trim(strn)
    call mksym(lc,slabl,strn,ips, iclasst,nclasst,npgrp,ngrp,oag,osymgr,oics,oistab,ngmx,faithful=.true.) !closed group with distinct rotations (gensym faithful in m_symfind)
    if(ipr10) write(stdo,"(a)") 'SpaceGroupSym of Lattice: ========end =========================== '
    allocate(symops,source=osymgr)
    allocate(ag,source=oag)
    AFmodeBlock: block
      real(8):: osymgrAll(3,3,ngmx),oagAll(3,ngmx)
      integer:: oicsAll(nbas),oistabAll(nbas,ngmx),ipsAF(nbas),iga,igall,nclassAll,ig
      AFmode=len_trim(symgaf)>0 
      if(AFmode) then
         strn2=trim(strn)//' '//trim(symgaf)
         if(ipr10) then
            write(stdo,*)
            write(stdo,"(a)") 'Add SpaceGroupSym ops by AF symmetry===start========= '
            write(stdo,"(a)") 'AF: Antiferro mode: SYMGRPAF    = '//trim(symgaf)
            write(stdo,"(a)") 'AF:                 SYMGRPAF all= '//trim(strn2)
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
         if(master_mpi) call pshpr(50)
         allocate(iclasstAll(nbas))
         if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice+AF: ========start========================== '
         call mksym(lc,slabl,strn2,ipsAF, iclasstAll,nclassAll, npgrpAll,ngrpTotal,oagAll,osymgrAll,oicsAll,oistabAll,ngmx) !Big symmetry for lattice+AF
         ngrpAF=ngrpTotal-ngrp
         if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice+AF: ========end========================== '
         if(master_mpi) call poppr()
!         do ig=1,ngrpTotal; write(stdo,ftox)'symall ig=',ig,ftof(reshape(osymgrAll(:,:,ig),[9]),2),' ',oagAll(:,ig);   enddo
         iga=ngrp
         do igall=1,ngrpTotal !Pick up symmetry by AF
            if(any( [( sum(abs(symops(:,:,ig)-osymgrALL(:,:,igall)))+sum(abs(ag(:,ig)-oagAll(:,igall)))<tol,ig=1,ngrp )]  )) cycle
            iga=iga+1
            symops(:,:,iga)= osymgrALL(:,:,igall)
            ag(:,iga)      = oagALL(:,igall)
            oistab(:,iga)  = oistabAll(:,igall)
         enddo
         if(iga/=ngrpAF+ngrp) call rxiii('ngrpAF+nggp/=ngrpTotal',ngrp,ngrpAF,iga)
      else
         ngrpTotal=ngrp
         ngrpAF=0
      endif
    endblock AFmodeBlock
    MiatTiatDlmm:block
      allocate(miat(nbas,ngrpTotal),tiat(3,nbas,ngrpTotal),invgx(ngrpTotal),shtvg(3,ngrpTotal),&
           dlmm(-lmxax:lmxax,-lmxax:lmxax,0:lmxax,ngrpTotal))
      ! The translations ag are given to mptauof (2026-10-02 05:42, step S4): shtvg = ag, and the inverse is found by rotation
      ! and translation. mptauof used to search a translation of its own, the same one for two operations with the same rotation.
      call            mptauof(symops,             ngrp,  plat,nbas,rv_a_opos,iclasst,   miat,tiat,invgx,shtvg, &
           ag=ag(:,1:ngrp))  !for ig=1,ngrp
      if(AFmode) call mptauof(symops(:,:,ngrp+1:),ngrpAF,plat,nbas,rv_a_opos,iclasstAll, & !  ig=ngrp+1,ngrpAF
           miat(:,ngrp+1:),tiat(:,:,ngrp+1:),invgx(ngrp+1:),shtvg(:,ngrp+1:),afmode, ag=ag(:,ngrp+1:ngrp+ngrpAF)) !mapping of sites by spacegrope ops
!      write(stdo,ftox)'mmmmm iclasst=',iclasst
!      do ig=1,ngrp
!         write(stdo,ftox)'mmmm ig=',ig, 'miat=',miat(1:nbas,ig)
!      enddo
!      write(stdo,ftox)'mmmmm iclasstAll=',iclasstAll
!      do ig=ngrp+1,ngrp+ngrpAF
!         write(stdo,ftox)'mmmm ig=',ig,'symops=',ftof(reshape(symops(:,:,ig),[9])),'miat=',miat(1:nbas,ig)
!      enddo
      call rotdlmm(symops,ngrpTotal, lmxax+1, dlmm) ! Get rotation matrix Dlmm in real spherical harmonics.  !for sigm mode, dlmm needed.
      if(ipr10) write(stdo,*)
    endblock MiatTiatDlmm
    call tcx('m_mksym_init')
  end subroutine m_mksym_init
  subroutine mksym(modeAddinversion,slabl,ssymgr,iv_a_oips, iclass,nclass,npgrp,nsgrp,rv_a_oag,rv_a_osymgr,iv_a_oics,iv_a_oistab,ngmxs,faithful)! Setup symmetry group. Split species into classes, Also assign class labels to each class
    use m_lmfinit,only: nbas,nspec,alat=>lat_alat,symgaf
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat,rv_a_opos
    use m_symfind,only: gensym,symfind_json,ngmx,ngnmx
    use m_ext,only: sname
    use m_symderive,only: grpgen,splcls,symtbl
    use m_symop_util,only: asymop
    use m_lgunit,only: stdo
    use m_ftox
    use m_mpi,only: master_mpi
    implicit none
    intent(in)::   modeAddinversion,slabl,ssymgr,iv_a_oips
    intent(out)::                                            iclass,nclass,npgrp,nsgrp,rv_a_oag,rv_a_osymgr,iv_a_oics,iv_a_oistab
    !i modeAddinversion  : 
    !i           =0  Not add inversion
    !i           =1  Add inversion to point group. Make additionally ag,istab for extra operations, using -g for rotation part; see Remarks
    !i slabl : species labels
    !i ssymgr: string containing symmetry group generators.
    !i           if ssymgr contains 'find', mksym will add basis atoms as
    !i           needed to guarantee generators are valid, and generate
    !i           internally any additonal group operations needed to
    !i           complete the space group.
    !r Remarks
    !r   In certain cases the inversion operation may be added to the space group, for purposes of k integration.  This is permissible when the
    !r   hamiltonian has the form h(-k) = h*(k).  In that case, the eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
    !r
    !r   Also, the Green's functions are related G(-k) = Gtranspose(k). Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
    !r   then G(-k) = Gtranspose(k), and the same (g,ag) information is needed for either rotation.
    integer :: modeAddinversion,nsgrp,npgrp,ibas,iwdummy1(1),idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nclass,ngen, nggen,incli
    integer,intent(in):: ngmxs ! size of the operation arrays (gensym fills at most ngmx=48 of them)
    integer:: iv_a_oips(nbas),iclass(nbas),ifind, iv_a_oics(nbas),iv_a_oistab(ngmxs*nbas)
    character(8) :: slabl(*),ssymgr*(*)
    character(1000) :: gens
    real(8) :: gen(3,3,ngnmx), rv_a_oag(3,ngmxs),rv_a_osymgr(3,3,ngmxs)
    integer,allocatable ::  iv_a_onrc (:), iv_a_oipc(:) 
    logical:: symfind
    logical,optional,intent(in):: faithful
    logical:: lfaithful,usejson
    character(600):: fjson
    character(200):: why
    character(16):: envv
    character(80):: sg
    integer:: ib
    lfaithful = .false.
    if(present(faithful)) lfaithful = faithful
    ifind = index(ssymgr,'find')
    gens = ssymgr
    symfind = ifind>0
    if(ifind>0) gens= ssymgr(1:ifind-1)//' '//ssymgr(ifind+4:)
    if(master_mpi) write(stdo,*)' Generators except find: ',trim(gens)
    if(master_mpi) write(stdo,*)' Generators find or not: ',symfind
    ! Backend of the finder (2026-10-02 05:35, step S3 of MD/symmetry_spglib.md): the operations of symmetry.<sname>.json
    ! (symfind.py, spglib) when the file is there, SYMGRP is 'find' (the default) and this is the group of the crystal (not the
    ! lattice+AF group, faithful=.false.); otherwise, or when the file cannot be used yet (AF), gensym. Pure translations are
    ! used as operations from step S4 on (2026-10-02 05:46).
    ! Not with SYMGRPAF yet (step S5): gensym returns its generators in ssymgr, and the AF call of m_mksym_init builds the
    ! lattice+AF group from them; with 'find' left there it would take another path.
    ! ECALJ_SYMFIND=ecalj in the environment forces gensym (for comparisons).
    usejson = .false.
    if(lfaithful .and. trim(adjustl(ssymgr))=='find' .and. len_trim(symgaf)==0) then
       fjson = 'symmetry.'//trim(sname)//'.json'
       inquire(file=trim(fjson),exist=usejson)
       envv = ' '
       call get_environment_variable('ECALJ_SYMFIND',envv)
       if(trim(envv)=='ecalj') usejson=.false.
    endif
    if(usejson) then
       call symfind_json(fjson,nbas,[(slabl(iv_a_oips(ib)),ib=1,nbas)],rv_a_opos(:,1:nbas),plat,qlat,alat, ngmxs, &
            nsgrp,rv_a_osymgr,rv_a_oag,usejson,why)
       if(master_mpi.and..not.usejson) write(stdo,"(a)")' mksym: '//trim(fjson)//' not used: '//trim(why)//'; gensym'
    endif
    if(usejson) then
       ngen  = 0
       nggen = nsgrp
       call symtbl(0,nbas,rv_a_opos,rv_a_osymgr,rv_a_oag,nsgrp,qlat,iv_a_oistab) !site ib goes to istab(ib,ig), as gensym returns it
       if(master_mpi) then
          write(stdo,"(a,i0,a)")' mksym: space group from '//trim(fjson)//' (spglib), ',nsgrp,' operations'
          write(stdo,"(' symfind_json: ig group ops (:vector means translation in cartesian)')")
          do ig = 1, nsgrp
             call asymop(rv_a_osymgr(:,:,ig),rv_a_oag(1,ig),':',sg)
             write(stdo,'(i5,2x,a)') ig,trim(sg)
          enddo
          write(stdo,"(a)")' symfind_json: site permutation table for group operations ...'
          write(stdo,"('  ib/ig:',48i3)")  [(ig,ig=1,nsgrp)]
          do ib = 1, nbas
             write(stdo,"(i7,':',48i3)") ib,(iv_a_oistab(ib+nbas*(ig-1)), ig=1,nsgrp)
          enddo
       endif
    else
       call gensym(slabl,gens,symfind,nbas,nspec,ngmx,plat,plat,rv_a_opos(:,1:nbas),iv_a_oips, & !Generate space group ops
            nsgrp, rv_a_osymgr,rv_a_oag, ngen,gen,ssymgr, nggen,isym,iv_a_oistab,lfaithful)
    endif
    if(nggen>ngmxs) call rx('mksym: nggen>ngmxs')
    incli = -1
    npgrp = nsgrp
    if(modeAddinversion /= 0) then !Add inversion to point group
       ngen = ngen+1
       gen(:,:,ngen) = reshape([-1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0],[3,3])
       call pshpr(iprint()-40)
       call grpgen(gen(1,1,ngen),1, rv_a_osymgr,npgrp, ngmxs)
       call poppr
       incli = npgrp-nsgrp
    endif  ! Printout of symmetry operations !    if(master_mpi) write(stdo,ftox)'  mksym: found ',nsgrp,' space group operations'
    if(master_mpi.and.nsgrp/=npgrp) write(stdo,ftox) &
         '    adding inversion gives',npgrp,' operations for generating k points; enforce real for dmatu for LDA+U'
    if(master_mpi.and.incli == -1) write(stdo,*)'  no attempt to add inversion symmetry'
    allocate(iv_a_onrc(nspec))
    allocate(iv_a_oipc,source=iv_a_oips(1:nbas))
    call splcls(rv_a_opos,nbas,nsgrp,iv_a_oistab,nspec,slabl,nclass,iv_a_oipc,iv_a_oics,iv_a_onrc) !Split species into classes
    !                                                   ibas ==> iclass=ipc(ibas) ==> ispec=ics(iclass)
    !  allocate(iv_a_oistab(nsgrp*nbas))
    call symtbl(1, nbas, rv_a_opos , rv_a_osymgr, rv_a_oag, nsgrp, qlat, iv_a_oistab)
    iclass(1:nbas)=iv_a_oipc(1:nbas) 
  end subroutine mksym
end module m_mksym
