module m_mksym
  use m_mksym_util,only: gensym,grpgen,symtbl
  public :: m_mksym_init, &
       ! ibas ==> iclass=ipc(ibas),iv_a_oipc ==> ispec=ics(iclass),iv_a_oics
       rv_a_oag,   iv_a_oics, iv_a_oipc , lat_npgrp, lat_nsgrp, &
       rv_a_osymgr,iv_a_oistab, ctrl_nclass, iclasstaf_, symops_af_,ag_af_,ngrpaf_

  integer, allocatable,protected ::  iv_a_oics (:)
  integer,  allocatable,protected ::  iv_a_oipc(:)
  real(8) , allocatable,protected ::  rv_a_oag (:)
  real(8) , allocatable,protected ::  rv_a_osymgr (:)
  integer , allocatable,protected ::  iv_a_oistab (:)
  integer,allocatable,protected:: iclasstaf_(:)
  real(8),allocatable,protected:: symops_af_(:,:,:), ag_af_(:,:)
  integer,protected:: lat_npgrp, lat_nsgrp, ctrl_nclass, ngrpaf_

  private
  integer, allocatable,protected ::  iv_a_onrc (:)
contains
  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine m_mksym_init(prgnam) !latsym(prgnam)
    use m_lmfinit,only: v_ssite,nbas,sstrnsymg,ctrl_noinv, &
         symgaf,iv_a_oips,slabl,mxspec,procid,master
    ! v_ssite%pos modified
    use m_lattic,only: rv_a_opos,m_lattic_init
    !! Driver for calling mksymaf and mksym
    !-------------
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
    !        ATOM=O POS=  .5   .5   .5
    !        ATOM=O POS= 1.5  1.5  1.5
    !------------------------------------------------------------
    implicit none
    character,intent(in)::  prgnam*(*)
    integer:: ibas,lc,j,iprint
    integer,parameter::recln=511
    logical ::cmdopt,ipr10=.false.
    character strn*(recln),strn2*(recln),outs(recln)
    call tcn('m_mksym_init')
    !! lattice set up
    do ibas=1,nbas
       v_ssite(ibas)%pos = rv_a_opos(:,ibas)
    enddo
    !! mksym mksymaf generate symmetry operations; split species into classes  ---
    !! mksym may correct atomic position slightly (not tested well).
    strn = 'find'
    if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
    if (cmdopt('--nosym',7,0,outs) .OR. cmdopt('--pdos',6,0,outs)) strn = ' '
    !! when lmfgw, 1st digit of lc is zero--> no inversion added in mksym.
    lc = 10
    if ( .NOT. prgnam == 'LMFGWD') lc = 12 ! add inversion for (not lmfgw)
    if (prgnam=='LMFA' .OR. prgnam=='LMCHK') then
       lc=20
       if(ctrl_noinv /=1 ) lc = lc+1  !+1 means add inversion
    endif
    if(len_trim(symgaf)>0) lc=12 ! inversion allowed for AF case. Good?
    if( .NOT. (prgnam=='LMFA' .OR. prgnam=='LMCHK')) ipr10= iprint()>10 !this is only for master
    if(len_trim(symgaf)>0) then
       if(ipr10) then
          write(6,*)
          write(6,"(a)")       ' AF: ======================================== '
          write(6,"(a)")       ' AF: Antiferro mode: SPGGRAF='//trim(symgaf)
          write(6,"(a)")       ' AF:  (neglct waring in GENSYM) '
          do j=1,nbas
             write(6,"(a,2i3)") ' AF:  ibas,AF=',j,v_ssite(j)%iantiferro
          enddo
       endif
       strn2=trim(strn)//' '//trim(symgaf)
       call mksymaf(v_ssite,iv_a_oips,nbas,procid==master,strn2,lc,slabl,mxspec)
       if(ipr10) write(6,"(a)") ' AF: ===== end of AF section================= '
       if(ipr10) write(6,"(a)")
    endif
    if(procid==master) call pshpr(60)
    call mksym(lc,slabl,strn,v_ssite,iv_a_oips)
    if(procid==master) call poppr
    call tcx('m_mksym_init')
  end subroutine m_mksym_init

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine mksymaf(v_ssite,iv_a_oips_in,nbas,imaster,strn2,lc,slabl,mxspec)!all input
    !     mksymaf sets data above used in sugw.F. a little complicated...
    use m_struc_def,only: s_site
    intent(in)::       v_ssite,iv_a_oips_in,nbas,imaster,strn2,lc,slabl,mxspec
    logical::                         imaster
    character strn2*(*)
    character(8) :: slabl(mxspec)
    type(s_site):: v_ssite(nbas)
    type(s_site),allocatable:: v_ssite2(:) !pos is rewritten in lmfp
    integer:: inumaf,iv_a_oips_in(nbas),j,k,nbas,  lc,       mxspec,ib
    integer,allocatable::iv_a_oips(:)
    allocate(v_ssite2(nbas),iv_a_oips(nbas))
    v_ssite2 = v_ssite
    iv_a_oips=iv_a_oips_in
    inumaf = 0
    do j=1,nbas
       do k=j,nbas
          if( v_ssite2(j)%iantiferro+v_ssite2(k)%iantiferro==0) then
             iv_a_oips(k) = iv_a_oips(j) !to drive mksymx
             inumaf=inumaf+1
             exit
          endif
       enddo
    enddo
    if(imaster) call pshpr(60)
    call mksym(lc,slabl,strn2,v_ssite2,iv_a_oips) !strn2 and v_ssite2 are used.
    if(imaster) call poppr()
    if(imaster) write(6,"(a)")' AF: mksym, generator= SYMGRP+SYMGRPAF= '//trim(strn2)
    ngrpaf_     = lat_nsgrp
    allocate(iclasstaf_(nbas),symops_af_(3,3,ngrpaf_),ag_af_(3,ngrpaf_))
    call dcopy ( ngrpaf_ * 9 , rv_a_osymgr , 1 , symops_af_ , 1 )
    call dcopy ( ngrpaf_ * 3 , rv_a_oag ,    1 , ag_af_ , 1 )
    if(imaster) write(6,"(a,i3)") ' AF: ngrpaf=',ngrpaf_
    do ib=1,nbas
       iclasstaf_(ib)=v_ssite2(ib)%class
    enddo
    deallocate(v_ssite2)
  end subroutine mksymaf

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs
  subroutine mksym(mode,slabl,ssymgr,ssite,iv_a_oips)
    use m_ftox
    use m_struc_def
    use m_lattic,only: lat_plat,rv_a_opos!,lat_dist
    use m_lmfinit,only: ctrl_nbas,ctrl_nspec,stdo!,lat_ldist
    !- Setup for symmetry group
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  : 1s digit
    !i           0  make space group only
    !i           1  add inversion to point group (see aginv below)
    !i           2  Same as 1, but make additionally ag,istab for extra
    !i              operations, using -g for rotation part; see Remarks
    !i           10s digit
    !i           0  do nothing about classes
    !i           1  split species into classes
    !i           2  Also assign class labels to each class
    !i           4  Assume class and species are the same.
    !i   slabl : species labels
    !i   ssymgr: string containing symmetry group generators.
    !i           if ssymgr contains 'find', mksym will add basis atoms as
    !i           needed to guarantee generators are valid, and generate
    !i           internally any additonal group operations needed to
    !i           complete the space group.
    !i
    ! o Inputs/outputs
    !i   sctrl :struct for program flow parameters; see routine uctrl
    !i     Elts read: nbas nspec lpgf
    !i     Stored:    nclass
    !i     Passed to: *
    ! o  slat  :struct containing information about the lattice
    ! o    Elts read: plat ldist dist
    ! o    Stored:    npgrp nsgrp oistab oag osymgr
    ! o          :nsgrp= number of space group operations.
    ! o          :npgrp= number of point group operations.  It may differ
    ! o                  from nsgrp if inversion is added.
    ! o          :istab= table of site permutations for each group op
    ! o          :ag   = translation part the group ops
    ! o          :symgr= point group (rotation) part of each group op
    ! o  ssite :struct containing site-specific information
    ! o    Elts read: *
    !i     Stored:    class clabel
    ! o  sarray:structure containing offsets to various arrays
    ! o    Elts read: opos oips npadl npadr
    ! o    Stored:    onrcp oipcp oics oipc oclabl onrc
    !o           oipc,oipcp,onrc,oics and oclabl, and
    !o           oipc:  pointer to class table, pad equiv to double pad
    !o           oipcp: pointer to class table, pad nequiv to double pad
    !o Outputs
    !o  ssymgr:  may be enlarged if it contains 'find' initially
    !o  slat     the following are set: npgrp nsgrp oistab oag osymgr
    !r Remarks
    !r   In certain cases the inversion operation may be added to the space
    !r   group, for purposes of k integration.  This is permissible when the
    !r   hamiltonian has the form h(-k) = h*(k).  In that case, the
    !r   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
    !r
    !r   Also, the Green's functions are related G(-k) = Gtranspose(k).
    !r   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
    !r   then G(-k) = Gtranspose(k), and the same (g,ag) information is
    !r   needed for either rotation.
    !u Updates
    !u   04 Jan 06 Enabled keyword 'point' in ssymgr
    !u             Generates extra group ops when generators create more
    !u             than the maximum expected
    !u   07 Oct 05 Put altered positions into site structure
    !u   03 Nov 01 Adapted to altered gensym
    !u   26 Jan 01 Added new '2' option to 1s digit of mode
    !u   31 May 00 revised meaning of mode
    ! ----------------------------------------------------------------------
    implicit none
    type(s_site)::ssite(*)
    character(8) :: slabl(*),ssymgr*(*)
    integer :: mode,nsgrp,npgrp,ibas
    integer ::iwdummy ,iwdummy1(1)
    logical :: T,F,cmdopt0,ltmp
    integer:: idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nbas, &
         nbas0,nclass,ngen,ngnmx,nspec,usegen, &
         nggen,ngmx,incli, oiwk , aginv  !ldist,
    integer,allocatable :: nrspc_iv(:)
    real(8) ,allocatable :: pos2_rv(:,:)
    integer ,allocatable :: ips2_iv(:)
    integer,allocatable:: iv_a_tmp(:)
    parameter (T=.true., F=.false., ngnmx=10)
    character(120) :: gens,strn(72)
    double precision :: gen(9,ngnmx),plat(3,3),qlat(3,3),xx !,dist(3,3)
    integer:: i_copy_size,i_data_size,i_spackv
    integer:: iv_a_oips(:)
    nbas =ctrl_nbas
    nspec=ctrl_nspec
    !      nsite=ctrl_nsite
    plat =lat_plat
    !      ldist=lat_ldist           ! ldist mode need to be recovered if necessary
    !      dist=lat_dist
    !      i_copy_size=size(lat_dist)
    !      call dcopy(i_copy_size,lat_dist,1,dist,1)
    !      stdo = globalvariables%stdo
    ngmx = 48
    ! ... Re-entry when ngmx was increased
5   continue
    ! ... For memory allocation of class data, assume initially nclass=nbas
    if(allocated(rv_a_oag)) deallocate(rv_a_oag,rv_a_osymgr,iv_a_oipc,iv_a_oics)
    allocate( rv_a_oag(3*ngmx)    )
    allocate( rv_a_osymgr(9*ngmx) )
    allocate( iv_a_oipc(nbas)  )
    allocate( iv_a_oics(nbas)  )
    allocate( nrspc_iv(nbas) )
    call words(ssymgr,ngen)
    j1 = 1
    idest = 1
    usegen = 2
    gens = ' '
    ltmp = .false.
    do  ig = 1, ngen
       call word(ssymgr,ig,j1,j2)
       if (ssymgr(j1:j2) == 'find') then
          usegen = 0
          !        else if (ssymgr(j1:j2) .eq. 'point') then
          !          ltmp = .true.
       else
          call strncp(gens,ssymgr,idest,j1,j2-j1+2)
          idest = idest+j2-j1+2
       endif
    enddo
    !      if (ltmp) then
    !        call strncp(gens,'point',idest,1,5)
    !      endif
    ! --- Generate space group ---
    nbas0 = nbas
    if (cmdopt0('--fixpos')) call Rx('fixpos is removed. outside of lmf')
    ! ... When generating the group the basis may become enlarged ...
    if(allocated(iv_a_oistab)) deallocate(iv_a_oistab) !lat_iv_a_oistab
    allocate(iv_a_oistab(abs((ngmx+1)*nbas)))
    allocate(ips2_iv(ngmx*nbas))
    allocate(pos2_rv(3,ngmx*nbas))
    ips2_iv(1:nbas)= iv_a_oips(1:nbas)
    call dcopy ( 3 * nbas , rv_a_opos , 1 , pos2_rv , 1 )
    call gensym ( slabl , gens , usegen , t , f , f , nbas &
         , nspec , ngmx , plat , plat , pos2_rv , ips2_iv& ! & , ldist , dist
         , nrspc_iv , nsgrp , rv_a_osymgr , rv_a_oag , ngen , gen , ssymgr &
         , nggen , isym , iv_a_oistab )
    if (nbas > nbas0) call rxs('gensym: the basis was enlarged.',' Check group operations.')
    if (nggen > nsgrp) then
       write(stdo,ftox)' MKSYM (warning): generators create more than ngmx=',ngmx,' group ops ...'
       ngmx = ngmx*16
       if (allocated(pos2_rv)) deallocate(pos2_rv)
       if (allocated(ips2_iv)) deallocate(ips2_iv)
       if (allocated(nrspc_iv)) deallocate(nrspc_iv)
       goto 5
    endif
    ! --- Add inversion to point group ---
    incli = -1
    npgrp = nsgrp
    if (mod(mode,10) /= 0) then
       ngen = ngen+1
       gen(:,ngen) = [-1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0]
       call pshpr(iprint()-40)
       call grpgen ( gen ( 1 , ngen ) , 1 , rv_a_osymgr , npgrp , ngmx  )
       call poppr
       incli = npgrp-nsgrp
    endif
    ! --- Printout of symmetry operations ---
    write(stdo,ftox)' MKSYM: found ',nsgrp,' space group operations'
    if(nsgrp/=npgrp) write(stdo,ftox)' MKSYM: adding inversion gives ',npgrp,' operations'
    if(incli == -1) write(stdo,*)'  no attempt to add inversion symmetry'
    if(mod(mode/10,10) == 0) goto 100
    ! Split species into classes : ibas ==> iclass=ipc(ibas) ==> ispec=ics(iclass)
    if(allocated(iv_a_onrc)) deallocate(iv_a_onrc)
    allocate(iv_a_onrc(abs(nspec)))
    iv_a_oipc=iv_a_oips(1:nbas)
    call splcls ( mod ( mode / 10 , 10 ) .eq.4 , rv_a_opos , nbas &
         , nsgrp , iv_a_oistab , nspec , slabl , nclass , iv_a_oipc , &
         iv_a_oics , iv_a_onrc )
    !   ... Reallocate arrays as permanent arrays
    i_data_size=size(iv_a_oics); allocate(iv_a_tmp(i_data_size))
    iv_a_tmp=iv_a_oics; deallocate(iv_a_oics)
    i_data_size=min(i_data_size,nclass); allocate(iv_a_oics(nclass))
    iv_a_oics(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)
    ! ... Remake istab
    if (allocated(iv_a_oistab)) deallocate(iv_a_oistab)
    allocate(iv_a_oistab(abs(nsgrp*nbas)))
    call dinv33(plat,1,qlat,xx)
    call symtbl(1, nbas, iwdummy1, rv_a_opos , rv_a_osymgr, rv_a_oag, nsgrp, qlat, iv_a_oistab)
    ! ... poke ipc into ssite
    do ibas=1,nbas
       ssite(ibas)%class = iv_a_oipc(ibas)
    enddo
    ! --- Cleanup: poke class and symmetry info into structures ---
    ctrl_nclass=nclass
100 continue
    lat_npgrp=npgrp
    lat_nsgrp=nsgrp
  end subroutine mksym
end module m_mksym