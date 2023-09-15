module m_mksym !crystal symmetry data given by call m_mksym_init.  nbas (atomic sites)-> nspec (species) -> nclass (class)
  use m_mpitk,only: master_mpi
  use m_lgunit,only: stdo
  public :: m_mksym_init, &
       rv_a_oag,   iv_a_oics, lat_npgrp, lat_nsgrp, &
       rv_a_osymgr,iv_a_oistab, ctrl_nclass,iclasst, &
       iclasstAF, symopsAF, agAF, ngrpAF !for antiferro symmetry
  !antiferro sym y = matmul(sympsAF(:,:,ig),x)+ ag_af(:,ig), ig=1,ngrpAF
  integer,  allocatable,protected :: iv_a_oics(:)     ! ispec= ics(iclass) gives spec for iclass.
  real(8) , allocatable,protected :: rv_a_osymgr (:,:,:) ! point operation for 1,... nsgrp  !(additioanl to npgrp)
  real(8) , allocatable,protected :: rv_a_oag (:)        ! translation for 1,...nsgrp
  integer , allocatable,protected :: iv_a_oistab (:)     ! j= istab(i,ig): site i is transformed into site j by grp op ig
  integer,allocatable,protected   :: iclasst(:)          !class information,
  integer,allocatable,protected:: iclasstAF(:) !=== AntiFerro class information 
  real(8),allocatable,protected:: symopsAF(:,:,:), agAF(:,:) 
  integer,protected::&
       lat_nsgrp,  & !number of space group symmetry
       ctrl_nclass,& !number of equivalence class
       lat_npgrp,  & !number of point group operation with adding iversion. Used in mkqp.f90
       ngrpAF       !antiferro symopsAF for (symopsAF,agAF)
  logical:: AFmode
contains
  subroutine m_mksym_init(prgnam) ! Driver for calling mksymaf and mksym
    use m_lmfinit,only: nspec,nbas,sstrnsymg,addinv, symgaf,iv_a_oips,slabl,mxspec,procid,master,iantiferro
    use m_lattic,only: rv_a_opos,m_lattic_init,rv_a_opos
    ! note: Antiferro symmetric self-energy mode. Feb2021 recovered
    !! To obtain self-consistency, it may be useful to keep AF condition during iteration.
    !! We need to set SYMGRPAF (still developing...)
    ! sample for nio
    !     SYMGRPAF i:(1,1,1) ! translation + inversion makes antiferro symmetry.
    !     SYMGRP r3d   !this keeps spin axis. (probably larger symmetry allowed for SO=0)
    ! TRUC   ALAT={a} PLAT= 0.5 0.5 1.0  0.5 1.0 0.5  1.0 0.5 0.5
    !        NBAS= 4  NSPEC=3
    ! ITE    ATOM=Niup POS=  .0   .0   .0   AF=1    <--- AF symmetric pair
    !        ATOM=Nidn POS= 1.0  1.0  1.0   AF=-1   <--- 
    !         ATOM=Mnup POS=  .0   .0   .0   AF=2   (if another AF pairs)
    !         ATOM=Mndn POS= 1.0  1.0  1.0   AF=-2 
    !        ATOM=O POS=  .5   .5   .5
    !        ATOM=O POS= 1.5  1.5  1.5
    !------------------------------------------------------------
    implicit none
    character,intent(in)::  prgnam*(*)
    integer:: ibas,lc,j,iprint
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
    AFsymblock :block
      integer:: inumaf,iv_a_oipsAF(nbas),j,k, mxspec,ib
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
         call mksym(lc,slabl,strn2,iv_a_oipsAF,iclasstAF) !lat_nsgrp returned !strn2 and v_ssite2 are used.
         ngrpAF    = lat_nsgrp !number of spg symmetry + af symmetry.
         allocate(symopsAF(3,3,ngrpAF),agAF(3,ngrpAF))
         if(master_mpi) call poppr()
         if(master_mpi) write(stdo,"(a)")'AF: mksym, generator= SYMGRP+SYMGRPAF= '//trim(strn2)
         call dcopy ( ngrpAF * 9 , rv_a_osymgr , 1 , symopsAF , 1 )
         call dcopy ( ngrpAF * 3 , rv_a_oag ,    1 , agAF , 1 )
         if(master_mpi) write(stdo,"(a,i3)")'AF: ngrpaf=',ngrpAF
         if(ipr10) write(stdo,"(a)") 'SpaceGroupSym of AF:      ========end =========================== '
         if(ipr10) write(stdo,"(a)")
      endif
    endblock AFsymblock
    if(procid==master) call pshpr(60)
    allocate(iclasst(nbas))
    if(ipr10) write(stdo,"(a)")      'SpaceGroupSym of Lattice: ========start========================== '
    call mksym(lc,slabl,strn,iv_a_oips,iclasst)
    if(ipr10) write(stdo,"(a)") 'SpaceGroupSym of Lattice: ========end =========================== '
    if(ipr10) write(stdo,*)
    if(procid==master) call poppr
    call tcx('m_mksym_init')
  end subroutine m_mksym_init
  
  subroutine mksym(mode,slabl,ssymgr,iv_a_oips,iclass)! Setup symmetry group.  Split species into classes, Also assign class labels to each class
    use m_lmfinit,only: nbas,stdo,nspec               ! outputs are allocated as module variables
    use m_lattic,only: lat_plat,rv_a_opos
    use m_mksym_util,only: gensym,grpgen,symtbl,splcls
    use m_ftox
    implicit none
    intent(in)::   mode,slabl,ssymgr,iv_a_oips
    intent(out)::                               iclass !class index
    !i Inputs
    !i   mode  : 
    !i           =0  Not add inversion
    !i           =1  Add inversion to point group 
    !i               Make additionally ag,istab for extra
    !i               operations, using -g for rotation part; see Remarks
    !i   slabl : species labels
    !i   ssymgr: string containing symmetry group generators.
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
    integer :: mode,nsgrp,npgrp,ibas,iwdummy1(1)
    integer:: idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nclass,ngen, nggen,ngmx,incli
    character(8) :: slabl(*),ssymgr*(*)
    logical :: cmdopt0
    real(8),allocatable :: pos2_rv(:,:)
    integer,allocatable :: ips2_iv(:)
    integer,parameter:: ngnmx=10
    character(120) :: gens
    real(8) :: gen(3,3,ngnmx),plat(3,3),qlat(3,3),xx
    integer:: iv_a_oips(:),iclass(nbas),ifind
    integer, allocatable ::  iv_a_onrc (:)
    integer, allocatable ::  iv_a_oipc(:) !class for lmaux                     maybe= iclasst
    logical:: symfind
    plat = lat_plat
    ngmx = 48
    if(allocated(rv_a_oag))  deallocate(rv_a_oag,rv_a_osymgr,iv_a_oics)
    if(allocated(iv_a_oipc)) deallocate(iv_a_oipc) !within subroutine
    allocate( rv_a_oag(3*ngmx)      )
    allocate( rv_a_osymgr(3,3,ngmx) )
    allocate( iv_a_oipc(nbas)  )
    allocate( iv_a_oics(nbas)  )
    ifind = index(ssymgr,'find')
    !write(stdo,*)'ifind=',ifind
    !write(stdo,*)'###'//ssymgr(ifind:ifind+3)//'###'
    gens = merge(ssymgr(1:ifind-1)//' '//ssymgr(ifind+4:),ssymgr, ifind>0)
    if(master_mpi) write(stdo,*)' Generators except find=',trim(gens)
    if(.not. allocated(iv_a_oistab)) allocate(iv_a_oistab(ngmx*nbas))
    symfind = ifind>0
    call gensym(slabl,gens,symfind,nbas,nspec,ngmx,plat,plat,rv_a_opos(:,1:nbas),iv_a_oips(1:nbas), & !Generate space group ops
         nsgrp,rv_a_osymgr,rv_a_oag, ngen,gen,ssymgr, nggen,isym,iv_a_oistab)
!    if(master_mpi) write(stdo,ftox)' mksym: ng ng ngen =',nsgrp,nggen,ngen
!    if(nggen>nsgrp.and.master_mpi) write(stdo,ftox)' MKSYM(warning): nggen=',nggen,'> nsgrp=',nsgrp
    if(nggen>ngmx) call rx('mksym: nggen>ngmx')
    incli = -1
    npgrp = nsgrp
    if (mode /= 0) then !Add inversion to point group
       ngen = ngen+1
       gen(:,:,ngen) = reshape([-1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0],[3,3])
       call pshpr(iprint()-40)
       call grpgen(gen(1,1,ngen),1, rv_a_osymgr,npgrp, ngmx)
       call poppr
       incli = npgrp-nsgrp
    endif
    ! Printout of symmetry operations !    if(master_mpi) write(stdo,ftox)'  mksym: found ',nsgrp,' space group operations'
    if(master_mpi.and.nsgrp/=npgrp) write(stdo,ftox) &
         '    adding inversion gives',npgrp,' operations for generating k points; enforce real for dmatu for LDA+U'
    if(master_mpi.and.incli == -1) write(stdo,*)'  no attempt to add inversion symmetry'
    if(allocated(iv_a_onrc)) deallocate(iv_a_onrc)
    allocate(iv_a_onrc(nspec))
    iv_a_oipc=iv_a_oips(1:nbas)
    call splcls(rv_a_opos,nbas,nsgrp,iv_a_oistab,nspec,slabl,nclass,iv_a_oipc,iv_a_oics,iv_a_onrc) !Split species into classes
    !                                                   ibas ==> iclass=ipc(ibas) ==> ispec=ics(iclass)
    if (allocated(iv_a_oistab)) deallocate(iv_a_oistab)
    allocate(iv_a_oistab(nsgrp*nbas))
    call dinv33(plat,1,qlat,xx)
    call symtbl(1, nbas, rv_a_opos , rv_a_osymgr, rv_a_oag, nsgrp, qlat, iv_a_oistab)
    iclass(1:nbas)=iv_a_oipc(1:nbas) 
    ctrl_nclass=nclass
    lat_npgrp=npgrp
    lat_nsgrp=nsgrp
  end subroutine mksym
end module m_mksym
